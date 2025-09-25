#
# The TreeSort application
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run qw(run);
use Cwd;
use Clone;

my $script = Bio::KBase::AppService::AppScript->new(\&process_treesort, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;


sub preflight
{
   # Declare and assign local variables for the parameters passed to the preflight function.
   my($app, $app_def, $raw_params, $params) = @_;

   print STDERR "Pre-flight TreeSort ", Dumper($params, $app);

   return {
      cpu => 2,
      memory => "64G",
      runtime => 18000,
      storage => 0,
   };
}

sub process_treesort
{
   # Declare and assign local variables for the parameters passed to the process_treesort function.
   my($app, $app_def, $raw_params, $params) = @_;

   # Uncomment to troubleshoot the parameters.
   warn Dumper($app_def, $raw_params, $params);

   my $token = $app->token();
   my $ws = $app->workspace();

   # Create a temp directory for intermediate and result files.
   my $cwd = File::Temp->newdir(CLEANUP => 1);

   # Create an "input" subdirectory for the input FASTA file, etc.
   my $input_dir = "$cwd/input";
   -d $input_dir or mkdir $input_dir or die "Cannot mkdir $input_dir: $!";

   # Create a "working" subdirectory for the input and intermediate file(s).
   my $work_dir = "$cwd/work";
   -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
   
   # TODO: Are these needed?
   my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
   my $dat = { data_api => $data_api };
   my $sstring = encode_json($dat);

   # Clone the input parameters.
   my $params_to_app = Clone::clone($params);

   # Encode the input parameters as a job description JSON file.
   my $job_desc = "$cwd/jobdesc.json";
   open(JOB_DESC, ">", $job_desc) or die "Cannot write $job_desc: $!";
   print JOB_DESC JSON::XS->new->pretty(1)->encode($params_to_app);
   close(JOB_DESC);

   my $parallel = $ENV{P3_ALLOCATED_CPU};

   # Run the Python script that runs TreeSort.
   my @cmd = ("run_treesort", "-i", $input_dir, "-j", $job_desc, "-w", $work_dir);
   my $ok = run(\@cmd);

   # Was the command successful?
   if (!$ok)
   {
      die "Command failed: @cmd\n";
      exit 1;
   }

   # Map file extensions to BV-BRC file types.
   my %suffix_map = (aln => 'aligned_dna_fasta',
                     csv => 'csv',
                     pdf => 'pdf',
                     tre => 'nwk',
                     tsv => 'tsv');

   my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

   if (opendir(my $dh, $work_dir))
   {
      while (my $p = readdir($dh))
      {
         next if $p =~ /^\./;

         # Don't copy the descriptor.csv file.
         next if $p eq "descriptor.csv";

         # Use the p3 utility to copy the files in the work directory to the user's workspace.
         my @cmd = ("p3-cp", "-r", "-f", @suffix_map, "$work_dir/$p", "ws:" . $app->result_folder);
         print "@cmd\n";

         my $ok = IPC::Run::run(\@cmd);
         if (!$ok)
         {
            warn "Error $? copying output with @cmd\n";
         }
      } 
      closedir($dh);
   }
   else
   {
      warn "Output directory $work_dir does not exist\n";
   }
}