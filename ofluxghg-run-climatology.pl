 #!/usr/bin/perl
 # script for running OceanFlux Greenhouse Gases climatology

 # v1 28/01/2013  ability to run 2010 fluxes (default to use Takahashi for CO2 data, but is ready to use SOCAT data including correct filename handling). Utils tested and working on the Nephalae Cloud, , jams@pml.ac.uk, Plymouth Marine Laboratory.

 # TO-DO: 
 # high priority 
 # 1. ability to set random noise values (rmse) in the config file
 # 2. check schmidt, kinematic viscosity and drag coefficient equations in the python against TS requirements
 
 
 # lower priority
 # 0. allow user defined output in netcdf file ? 


use strict;
use warnings;
use POSIX; 
use Getopt::Long;
use Sys::Hostname;
use File::Spec;
use File::Basename;
use File::Path;
use AppConfig;
#use File::HomeDir; # not installed on the CERSAT Cloud

# subroutine prototypes
sub usage;

# set program variables
my $prog_dir = File::Spec->rel2abs(dirname($0));
my $sw_version = "001";
my $prog = basename($0);
my $func = 'main';
my $hostname = hostname;
my $debug=1; # debug mode

 # locate the users home directory
 # preferred approach, needs File::HomeDir but isn't installed
# my $home_dir = File::HomeDir->my_home; 
 # alternative approach to determine home directory
my $home_dir = $ENV{HOME};
print "home dir: $home_dir\n";

 # variables to be read from the config
my ($SRC_HOME);
my ($WINDU10, $SSTSKIN, $SSTFND, $SIGMA0, $SIG_WV_HT, $ICE, $PRESSURE, $PCO2, $SALINITY, $RAIN, $PCO2_DATA_SELECTION, $SALINITY_DATA_SELECTION, $RAIN_DATA_SELECTION, $SSTGRAD, $BIOLOGY, $SUSPENDED_PARTICLES, $ATLANTIC_OCEAN, $PACIFIC_OCEAN, $SOUTHERN_OCEAN, $INDIAN_OCEAN, $LONGHURST);
my ($WINDU10_PROD_NAME, $SSTSKIN_PROD_NAME, $SSTFND_PROD_NAME, $SIGMA0_PROD_NAME, $SIG_WV_HT_PROD_NAME, $ICE_PROD_NAME, $PCO2_PROD_NAME, $PCO2_SST_PROD_NAME, $VCO2_PROD_NAME, $SALINITY_PROD_NAME, $RAIN_PROD_NAME);
my ($RANDOM_NOISE_WINDU10, $RANDOM_NOISE_SSTSKIN, $RANDOM_NOISE_SSTFND, $RANDOM_NOISE_PCO2, $BIAS_WINDU10, $BIAS_SSTSKIN, $BIAS_SSTFND, $BIAS_PCO2, $BIAS_K, $BIAS_WINDU10_VALUE, $BIAS_SSTSKIN_VALUE, $BIAS_SSTFND_VALUE, $BIAS_PCO2_VALUE, $BIAS_K_VALUE, $BIAS_K_PERCENT, $BIAS_K_BIOLOGY_VALUE, $BIAS_K_WIND_VALUE);
my ($K_HO2006, $K_NIGHTINGALE2000, $KT_OCEANFLUXGHG, $K_WANNINKHOF1992, $K_WANNINKHOF2013, $K_MCGILLIS2001, $K_HO1997, $KD_OCEANFLUXGHG_BACKSCATTER, $KD_OCEANFLUXGHG_WIND, $KB_OCEANFLUXGHG, $KT_OCEANFLUXGHG_KD_WIND, $K_GENERIC, $K_GENERIC_SC, $K_GENERIC_A0, $K_GENERIC_A1,$K_GENERIC_A2,$K_GENERIC_A3, $KD_WEIGHTING, $KB_WEIGHTING, $KB_ASYMMETRY);
my ($RAPID, $EQUILIBRIUM, $BULK, $SST_GRADIENTS, $USE_SSTSKIN, $USE_SSTFND, $SALINE_SKIN_VALUE);
my ($OUTPUT_DIR);


 # values fixed for climatology runs - needed for rain studies
my ($BIAS_SSTSKIN_DUE_RAIN, $BIAS_SSTSKIN_DUE_RAIN_VALUE, $BIAS_SSTSKIN_DUE_RAIN_INTENSITY, $BIAS_SSTSKIN_DUE_RAIN_WIND, $RAIN_WET_DEPOSITION, $K_RAIN_LINEAR_HO1997, $K_RAIN_NONLINEAR_H2012);
 # SST skin modification by rain
#$BIAS_SSTSKIN_DUE_RAIN = 0;
#$BIAS_SSTSKIN_DUE_RAIN_VALUE = 0.0;
#$BIAS_SSTSKIN_DUE_RAIN_INTENSITY = 0.0;
#$BIAS_SSTSKIN_DUE_RAIN_WIND = 0.0;
 # wet deposition by rain switch
$RAIN_WET_DEPOSITION = 0;
$K_RAIN_LINEAR_HO1997 = 0;
$K_RAIN_NONLINEAR_H2012 = 0;

 # grab options
my ($help, $config_file, $pco2_dir_override, $output_dir_override, $process_layers_off, $year_start, $year_end);

 # defaults
$year_start = 2010;
$year_end = 2010;

GetOptions(
   'help|h!'      =>\$help,
   'config|c=s'   =>\$config_file,
   'pco2_dir_override|p=s' =>\$pco2_dir_override,
   'output_dir_override|o=s'  =>\$output_dir_override,
   'process_layers_off|l'  =>\$process_layers_off,
   'year_start|s=i'  =>\$year_start,
   'year_end|e=i' =>\$year_end,
);


 # display the help information
if($help) {
   usage();
   exit(0);
}

 # check for the configuration file
if(!$config_file)
{
   usage();
   die "($prog, $func) NOTE: Missing --config <file>, exiting.";
}



 #setting up configuration object
my $config = AppConfig->new( {
    CREATE=>1,
  } );

 # define instances
$config->define('src_home', { ARGS  => "=s"});
$config->define('rapid', { ARGS  => "=s"});
$config->define('equilibrium', { ARGS  => "=s"});
$config->define('bulk', { ARGS  => "=s"});
$config->define('sst_gradients', { ARGS  => "=s"});
$config->define('use_sstskin', { ARGS  => "=s"});
$config->define('use_sstfnd', { ARGS  => "=s"});
$config->define('saline_skin_value', { ARGS  => "=f"});
$config->define('windu10', { ARGS  => "=s"});
$config->define('sstskin', { ARGS  => "=s"});
$config->define('sstfnd', { ARGS  => "=s"});
$config->define('sigma0', { ARGS  => "=s"});
$config->define('sig_wv_ht', { ARGS  => "=s"});
$config->define('ice', { ARGS  => "=s"});
$config->define('ice_prod_name', { ARGS  => "=s"});
$config->define('pressure', { ARGS  => "=s"});
$config->define('pco2', { ARGS  => "=s"});
$config->define('salinity', { ARGS  => "=s"});
$config->define('rain', { ARGS  => "=s"});
$config->define('pco2_data_selection', { ARGS  => "=s"});
$config->define('salinity_data_selection', { ARGS  => "=s"});
$config->define('rain_data_selection', { ARGS  => "=s"});
$config->define('windu10_prod_name', { ARGS  => "=s"});
$config->define('sstskin_prod_name', { ARGS  => "=s"});
$config->define('sstfnd_prod_name', { ARGS  => "=s"});
$config->define('sigma0_prod_name', { ARGS  => "=s"});
$config->define('sig_wv_ht_prod_name', { ARGS  => "=s"});
$config->define('pco2_prod_name', { ARGS  => "=s"});
$config->define('pco2_sst_prod_name', { ARGS  => "=s"});
$config->define('vco2_prod_name', { ARGS  => "=s"});
$config->define('salinity_prod_name', { ARGS  => "=s"});
$config->define('rain_prod_name', { ARGS  => "=s"});
$config->define('salinity', { ARGS  => "=s"});
$config->define('sstgrad', { ARGS  => "=s"});
$config->define('biology', { ARGS  => "=s"});
$config->define('atlantic_ocean', { ARGS  => "=s"});
$config->define('pacific_ocean', { ARGS  => "=s"});
$config->define('southern_ocean', { ARGS  => "=s"});
$config->define('indian_ocean', { ARGS  => "=s"});
$config->define('longhurst', { ARGS  => "=s"});
$config->define('suspended_particles', { ARGS  => "=s"});
$config->define('random_noise_windu10', { ARGS  => "=s"});
$config->define('random_noise_sstskin', { ARGS  => "=s"});
$config->define('random_noise_sstfnd', { ARGS  => "=s"});
$config->define('random_noise_pco2', { ARGS  => "=s"});
$config->define('bias_windu10', { ARGS  => "=s"});
$config->define('bias_sstskin', { ARGS  => "=s"});
$config->define('bias_sstfnd', { ARGS  => "=s"});
$config->define('bias_pco2', { ARGS  => "=s"});
$config->define('bias_k', { ARGS  => "=s"});
$config->define('bias_windu10_value', { ARGS  => "=f"});
$config->define('bias_sstskin_value', { ARGS  => "=f"});
$config->define('bias_sstfnd_value', { ARGS  => "=f"});
$config->define('bias_pco2_value', { ARGS  => "=f"});
$config->define('bias_k_value', { ARGS  => "=f"});
$config->define('bias_k_percent', { ARGS  => "=f"});
$config->define('bias_k_biology_value', { ARGS  => "=f"});
$config->define('bias_k_wind_value', { ARGS  => "=f"});
$config->define('k_Ho2006', { ARGS  => "=s"});
$config->define('k_Nightingale2000', { ARGS  => "=s"});
$config->define('kt_OceanFluxGHG', { ARGS  => "=s"});
$config->define('k_Wanninkhof1992', { ARGS  => "=s"});
$config->define('k_Wanninkhof2013', { ARGS  => "=s"});
$config->define('k_McGillis2001', { ARGS  => "=s"});
$config->define('k_Ho1997', { ARGS  => "=s"});
$config->define('kd_OceanFluxGHG_backscatter', { ARGS  => "=s"});
$config->define('kd_OceanFluxGHG_wind', { ARGS  => "=s"});
$config->define('kb_OceanFluxGHG', { ARGS  => "=s"});
$config->define('kt_OceanFluxGHG_kd_wind', { ARGS  => "=s"});
$config->define('k_generic', { ARGS  => "=s"});
$config->define('k_generic_sc', { ARGS  => "=f"});
$config->define('k_generic_a0', { ARGS  => "=f"});
$config->define('k_generic_a1', { ARGS  => "=f"});
$config->define('k_generic_a2', { ARGS  => "=f"});
$config->define('k_generic_a3', { ARGS  => "=f"});
$config->define('kd_weighting', { ARGS  => "=f"});
$config->define('kb_weighting', { ARGS  => "=f"});
$config->define('kb_asymmetry', { ARGS  => "=f"});
 # rain specific options
$config->define('bias_sstskin_due_rain', { ARGS  => "=s"});
$config->define('bias_sstskin_due_rain_value', { ARGS  => "=f"});
$config->define('bias_sstskin_due_rain_intensity', { ARGS  => "=f"});
$config->define('bias_sstskin_due_rain_wind', { ARGS  => "=f"});
$config->define('rain_wet_deposition', { ARGS  => "=s"});
$config->define('k_rain_linear_Ho1997', { ARGS  => "=s"});
$config->define('k_rain_nonlinear_H2012', { ARGS  => "=s"});

$config->define('output_dir', { ARGS  => "=s"});
 #read the file
$config->file($config_file);


 # grabbing entries
 # checking all entries exist in the configuration file

if (!($SRC_HOME  = $config->src_home)){
   die "($prog, $func) Failed to determine SRC_HOME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RAPID = $config->rapid)){
   die "($prog, $func) Failed to determine RAPID entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($EQUILIBRIUM = $config->equilibrium)){
   die "($prog, $func) Failed to determine EQUILIBRIUM entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BULK = $config->bulk)){
   die "($prog, $func) Failed to determine BULK entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SST_GRADIENTS = $config->sst_gradients)){
   die "($prog, $func) Failed to determine SST_GRADIENTS entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($USE_SSTSKIN = $config->use_sstskin)){
   die "($prog, $func) Failed to determine USE_SSTSKIN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($USE_SSTFND = $config->use_sstfnd)){
   die "($prog, $func) Failed to determine USE_SSTFND entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SALINE_SKIN_VALUE = $config->saline_skin_value)){
   die "($prog, $func) Failed to determine SALINE_SKIN_VALUE entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($WINDU10  = $config->windu10)){
   die "($prog, $func) Failed to determine WINDU10 entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SSTSKIN = $config->sstskin)){
   die "($prog, $func) Failed to determine SSTSKIN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SSTFND = $config->sstfnd)){
   die "($prog, $func) Failed to determine SSTFND entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SIGMA0 = $config->sigma0)){
   die "($prog, $func) Failed to determine SIGMA0 entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SIG_WV_HT = $config->sig_wv_ht)){
   die "($prog, $func) Failed to determine SIG_WV_HT entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($ICE = $config->ice)){
   die "($prog, $func) Failed to determine ICE entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($ICE_PROD_NAME = $config->ice_prod_name)){
   die "($prog, $func) Failed to determine ICE_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($WINDU10_PROD_NAME  = $config->windu10_prod_name)){
   die "($prog, $func) Failed to determine WINDU10_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SSTSKIN_PROD_NAME = $config->sstskin_prod_name)){
   die "($prog, $func) Failed to determine SSTSKIN_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SSTFND_PROD_NAME = $config->sstfnd_prod_name)){
   die "($prog, $func) Failed to determine SSTFND_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SIGMA0_PROD_NAME = $config->sigma0_prod_name)){
   die "($prog, $func) Failed to determine SIGMA0_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SIG_WV_HT_PROD_NAME = $config->sig_wv_ht_prod_name)){
   die "($prog, $func) Failed to determine SIG_WV_HT_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($PCO2_PROD_NAME = $config->pco2_prod_name)){
   die "($prog, $func) Failed to determine PCO2_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($PCO2_SST_PROD_NAME = $config->pco2_sst_prod_name)){
   die "($prog, $func) Failed to determine PCO2_SST_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($VCO2_PROD_NAME = $config->vco2_prod_name)){
   die "($prog, $func) Failed to determine VCO2_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SALINITY_PROD_NAME = $config->salinity_prod_name)){
   die "($prog, $func) Failed to determine SALINITY_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RAIN_PROD_NAME = $config->rain_prod_name)){
   die "($prog, $func) Failed to determine RAIN_PROD_NAME entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($PRESSURE = $config->pressure)){
   die "($prog, $func) Failed to determine PRESSURE entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($PCO2 = $config->pco2)){
   die "($prog, $func) Failed to determine PCO2 entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SALINITY = $config->salinity)){
   die "($prog, $func) Failed to determine SALINITY entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RAIN = $config->rain)){
   die "($prog, $func) Failed to determine RAIN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($PCO2_DATA_SELECTION = $config->pco2_data_selection)){
   die "($prog, $func) Failed to determine PCO2_DATA_SELECTION entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SALINITY_DATA_SELECTION = $config->salinity_data_selection)){
   die "($prog, $func) Failed to determine SALINITY_DATA_SELECTION entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RAIN_DATA_SELECTION = $config->rain_data_selection)){
   die "($prog, $func) Failed to determine RAIN_DATA_SELECTION entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SSTGRAD = $config->sstgrad)){
   die "($prog, $func) Failed to determine SSTGRAD entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIOLOGY = $config->biology)){
   die "($prog, $func) Failed to determine BIOLOGY entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($ATLANTIC_OCEAN = $config->atlantic_ocean)){
   die "($prog, $func) Failed to determine ATLANTIC_OCEAN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($PACIFIC_OCEAN = $config->pacific_ocean)){
   die "($prog, $func) Failed to determine PACIFIC_OCEAN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SOUTHERN_OCEAN = $config->southern_ocean)){
   die "($prog, $func) Failed to determine SOUTHERN_OCEAN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($INDIAN_OCEAN = $config->indian_ocean)){
   die "($prog, $func) Failed to determine INDIAN_OCEAN entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($LONGHURST = $config->longhurst)){
   die "($prog, $func) Failed to determine LONGHURST entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($SUSPENDED_PARTICLES = $config->suspended_particles)){
   die "($prog, $func) Failed to determine SUSPENDED_PARTICLES entry from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RANDOM_NOISE_WINDU10 = $config->random_noise_windu10)){
   die "($prog, $func) Failed to determine RANDOM_NOISE_WINDU10 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RANDOM_NOISE_SSTSKIN = $config->random_noise_sstskin)){
   die "($prog, $func) Failed to determine RANDOM_NOISE_SSTSKI value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RANDOM_NOISE_SSTFND = $config->random_noise_sstfnd)){
   die "($prog, $func) Failed to determine RANDOM_NOISE_SSTFND value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RANDOM_NOISE_PCO2 = $config->random_noise_pco2)){
   die "($prog, $func) Failed to determine RANDOM_NOISE_PCO2 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_WINDU10 = $config->bias_windu10)){
   die "($prog, $func) Failed to determine BIAS_WINDU10 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTSKIN = $config->bias_sstskin)){
   die "($prog, $func) Failed to determine BIAS_SSTSKIN value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTFND = $config->bias_sstfnd)){
   die "($prog, $func) Failed to determine BIAS_SSTFND value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_PCO2 = $config->bias_pco2)){
   die "($prog, $func) Failed to determine BIAS_PCO2 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_K = $config->bias_k)){
   die "($prog, $func) Failed to determine BIAS_K value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_WINDU10_VALUE = $config->bias_windu10_value)){
   die "($prog, $func) Failed to determine BIAS_WINDU10_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTSKIN_VALUE = $config->bias_sstskin_value)){
   die "($prog, $func) Failed to determine BIAS_SSTSKIN_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTFND_VALUE = $config->bias_sstfnd_value)){
   die "($prog, $func) Failed to determine BIAS_SSTFND_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_PCO2_VALUE = $config->bias_pco2_value)){
   die "($prog, $func) Failed to determine BIAS_PCO2_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_K_VALUE = $config->bias_k_value)){
   die "($prog, $func) Failed to determine BIAS_K_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_K_PERCENT = $config->bias_k_percent)){
   die "($prog, $func) Failed to determine BIAS_K_PERCENT value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_K_BIOLOGY_VALUE = $config->bias_k_biology_value)){
   die "($prog, $func) Failed to determine BIAS_K_BIOLOGY_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_K_WIND_VALUE = $config->bias_k_wind_value)){
   die "($prog, $func) Failed to determine BIAS_K_WIND_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_HO2006 = $config->k_Ho2006)){
   die "($prog, $func) Failed to determine K_HO2006 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_NIGHTINGALE2000 = $config->k_Nightingale2000)){
   die "($prog, $func) Failed to determine K_NIGHTINGALE2000 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KT_OCEANFLUXGHG = $config->kt_OceanFluxGHG)){
   die "($prog, $func) Failed to determine KT_OCEANFLUXGHG value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_WANNINKHOF1992 = $config->k_Wanninkhof1992)){
   die "($prog, $func) Failed to determine K_WANNINKHOF1992 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_WANNINKHOF2013 = $config->k_Wanninkhof2013)){
   die "($prog, $func) Failed to determine K_WANNINKHOF2013 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_MCGILLIS2001 = $config->k_McGillis2001)){
   die "($prog, $func) Failed to determine K_MCGILLIS2001 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_HO1997 = $config->k_Ho1997)){
   die "($prog, $func) Failed to determine K_HO1997 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KD_OCEANFLUXGHG_BACKSCATTER = $config->kd_OceanFluxGHG_backscatter)){
   die "($prog, $func) Failed to determine KD_OCEANFLUXGHG_BACKSCATTER value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KD_OCEANFLUXGHG_WIND = $config->kd_OceanFluxGHG_wind)){
   die "($prog, $func) Failed to determine KD_OCEANFLUXGHG value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KB_OCEANFLUXGHG = $config->kb_OceanFluxGHG)){
   die "($prog, $func) Failed to determine KB_OCEANFLUXGHG value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KT_OCEANFLUXGHG_KD_WIND = $config->kt_OceanFluxGHG_kd_wind)){
   die "($prog, $func) Failed to determine KT_OCEANFLUXGHG_KD_WIND value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_GENERIC = $config->k_generic)){
   die "($prog, $func) Failed to determine K_GENERIC value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_GENERIC_SC = $config->k_generic_sc)){
   die "($prog, $func) Failed to determine K_GENERIC_SC value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_GENERIC_A0 = $config->k_generic_a0)){
   die "($prog, $func) Failed to determine K_GENERIC_A0 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_GENERIC_A1 = $config->k_generic_a1)){
   die "($prog, $func) Failed to determine K_GENERIC_A1 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_GENERIC_A2 = $config->k_generic_a2)){
   die "($prog, $func) Failed to determine K_GENERIC_A2 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_GENERIC_A3 = $config->k_generic_a3)){
   die "($prog, $func) Failed to determine K_GENERIC_A3 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KD_WEIGHTING = $config->kd_weighting)){
   die "($prog, $func) Failed to determine KD_WEIGHTING value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KB_WEIGHTING = $config->kb_weighting)){
   die "($prog, $func) Failed to determine KB_WEIGHTING value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($KB_ASYMMETRY = $config->kb_asymmetry)){
   die "($prog, $func) Failed to determine KB_ASYMMETRY value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTSKIN_DUE_RAIN = $config->bias_sstskin_due_rain)){
die "($prog, $func) Failed to determine BIAS_SSTSKIN_DUE_RAIN value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTSKIN_DUE_RAIN_VALUE = $config->bias_sstskin_due_rain_value)){
die "($prog, $func) Failed to determine BIAS_SSTSKIN_DUE_RAIN_VALUE value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTSKIN_DUE_RAIN_INTENSITY = $config->bias_sstskin_due_rain_intensity)){
die "($prog, $func) Failed to determine BIAS_SSTSKIN_DUE_RAIN_INTENSITY value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($BIAS_SSTSKIN_DUE_RAIN_WIND = $config->bias_sstskin_due_rain_wind)){
die "($prog, $func) Failed to determine BIAS_SSTSKIN_DUE_RAIN_WIND value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($RAIN_WET_DEPOSITION = $config->rain_wet_deposition)){
die "($prog, $func) Failed to determine RAIN_WET_DEPOSITION value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_RAIN_LINEAR_HO1997 = $config->k_rain_linear_Ho1997)){
die "($prog, $func) Failed to determine K_RAIN_LINEAR_HO1997 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($K_RAIN_NONLINEAR_H2012 = $config->k_rain_nonlinear_H2012)){
die "($prog, $func) Failed to determine K_RAIN_NONLINEAR_H2012 value from file:$config_file, not a valid configuration file so exiting..";
}
if (!($OUTPUT_DIR = $config->output_dir)){
   die "($prog, $func) Failed to determine the output_directory from file:$config_file, not a valid configuration file so exiting..";
}


 # adding in home directory to each path
$WINDU10 = $home_dir."/".$WINDU10;
$SSTSKIN = $home_dir."/".$SSTSKIN;
$SSTFND = $home_dir."/".$SSTFND;
$SIGMA0 = $home_dir."/".$SIGMA0;
$SIG_WV_HT = $home_dir."/".$SIG_WV_HT;
$PRESSURE = $home_dir."/".$PRESSURE;
$ICE = $home_dir."/".$ICE;


#########################################
 ## commenting out to test new SOCAT data
#$PCO2 = $home_dir."/".$PCO2;
#########################################


$RAIN = $home_dir."/".$RAIN;
$SALINITY = $home_dir."/".$SALINITY;
$SSTGRAD = $home_dir."/".$SSTGRAD;
#$BIOLOGY = $home_dir."/".$BIOLOGY; # comment out for CCI data
$OUTPUT_DIR = $home_dir."/".$OUTPUT_DIR;
#$SUSPENDED_PARTICLES = $home_dir."/".$SUSPENDED_PARTICLES; # comment out for CCI data
$ATLANTIC_OCEAN = $home_dir."/".$ATLANTIC_OCEAN;
$PACIFIC_OCEAN = $home_dir."/".$PACIFIC_OCEAN;
$SOUTHERN_OCEAN = $home_dir."/".$SOUTHERN_OCEAN;
$INDIAN_OCEAN = $home_dir."/".$INDIAN_OCEAN;

$LONGHURST = $home_dir."/".$LONGHURST;


 # detection of optional arguments from command line
 
 # detection of optional filelist for pco2 data
if ($pco2_dir_override){
   print STDOUT "\n($prog, $func) Using optional override for pCO2w data directory, \$PCO2 will be set to $pco2_dir_override (ie overriding both directory of pco2 data and selection in the configuration file).";
    # set to using SOCAT data
   #$PCO2_DATA_SELECTION = "socat_pco2";
   $PCO2 = $pco2_dir_override;
}
 # detection of optional output_dir override
if ($output_dir_override) {
   print STDOUT "\n($prog, $func) Using optional override for output directory, output_dir will be set to $output_dir_override/<year>/<month>/ (overriding the configuration file value).";
   $OUTPUT_DIR = $output_dir_override;
}

if ($process_layers_off) {
   print STDOUT "\n($prog, $func) Switching off generation of processing indicator layers. This reduces the processing time by apx 50% (switch: $process_layers_off).";
   $process_layers_off = 1;
} else {
   $process_layers_off = 0;
}


#print "contents of configuration: $config_file_string_metadata\n";

 # checks that only one k parameterisation is chosen, only some combinations are listed here
if ( (($K_HO2006 eq 'yes') and ($K_NIGHTINGALE2000 eq 'yes')) or (($K_HO2006 eq 'yes') and ($KT_OCEANFLUXGHG eq 'yes') and ($K_HO1997 eq 'yes')) or (($K_HO2006 eq 'yes') and ($K_WANNINKHOF1992 eq 'yes')) or (($K_HO2006 eq 'yes') and ($K_MCGILLIS2001 eq 'yes') and ($K_HO1997 eq 'yes')) or (($KT_OCEANFLUXGHG eq 'yes') and ($K_NIGHTINGALE2000 eq 'yes')) or (($KT_OCEANFLUXGHG eq 'yes') and ($K_WANNINKHOF1992 eq 'yes') and ($K_HO1997 eq 'yes')) or (($KT_OCEANFLUXGHG eq 'yes') and ($K_MCGILLIS2001 eq 'yes') and ($K_HO1997 eq 'yes')) or (($K_WANNINKHOF1992 eq 'yes') and ($K_MCGILLIS2001 eq 'yes')) or (($K_NIGHTINGALE2000 eq 'yes') and ($K_MCGILLIS2001 eq 'yes') and ($K_HO1997 eq 'yes')) or (($K_NIGHTINGALE2000 eq 'yes') and ($K_WANNINKHOF1992 eq 'yes')) ){
  die "($prog, $func) Only one k parameterisation choice is possible, file:$config_file is not a valid configuration file so exiting..";
}

 # need to actually choose a k parameterisation
if ( ($K_HO2006 eq 'no') and ($K_NIGHTINGALE2000 eq 'no') and ($KT_OCEANFLUXGHG eq 'no') and ($K_WANNINKHOF1992 eq 'no') and ($K_WANNINKHOF2013 eq 'no') and ($K_MCGILLIS2001 eq 'no') and ($K_HO1997 eq 'no') and ($KD_OCEANFLUXGHG_BACKSCATTER eq 'no') and ($KB_OCEANFLUXGHG eq 'no') and ($K_GENERIC eq 'no') and ($KD_OCEANFLUXGHG_WIND eq 'no') and ($KT_OCEANFLUXGHG_KD_WIND eq 'no') and ($RAIN_WET_DEPOSITION eq 'no')){
   die "($prog, $func) No k parameterisation has been chosen, file:$config_file is not a valid configuration file so exiting..";
}

if ($K_GENERIC eq 'yes'){
   if (($K_GENERIC_SC != 600) and ($K_GENERIC_SC != 660)){
      die "($prog, $func) The 'k_generic' option has been selected.  Only schmidt numbers of 600 and 660 are allowed, file:$config_file is not a valid configuration file so exiting..";
   }
}

 # checking Rapid/equilibrium model selection
my $flux_model=0;
if ( (($RAPID eq 'yes') and ($EQUILIBRIUM eq 'yes')) or (($RAPID eq 'yes') and ($BULK eq 'yes')) or (($EQUILIBRIUM eq 'yes') and ($BULK eq 'yes'))) {
   die "($prog, $func) More than one flux model (RAPID, EQUILIBRIUM or BULK) selected in configuration file ($config_file), only one can be selected so this is not a valid configuration file so exiting..";
}
if ($RAPID eq 'yes'){ $flux_model = 1;}
elsif ($EQUILIBRIUM eq 'yes'){ $flux_model = 2;}
elsif ($BULK eq 'yes'){ $flux_model = 3;}
else {
   die "($prog, $func) Neither RAPID, EQUILIBRIUM or BULK models have been selected in configuration file ($config_file), not a valid configuration file so exiting..";
}

if ($debug ==1){
   print "\n\n($prog, $func) Flux model: $flux_model\n\n";
}

 # checking sst_gradients selection
my $sst_gradients= -1;
if ($SST_GRADIENTS eq 'yes'){ $sst_gradients = 1;}
elsif ($SST_GRADIENTS eq 'no'){ $sst_gradients = 0;}
else {
   die "($prog, $func) SST_GRADIENTS selection needs to be yes or no in configuration file ($config_file), not a valid configuration file so exiting..";
}

 # checking use_sstskin selection
my $use_sstskin= -1;
if ($USE_SSTSKIN eq 'yes'){ $use_sstskin = 1;}
elsif ($USE_SSTSKIN eq 'no'){ $use_sstskin = 0;}
else {
   die "($prog, $func) USE_SSTSKIN selection needs to be yes or no in configuration file ($config_file), not a valid configuration file so exiting..";
}

 # checking use_sstfnd selection
my $use_sstfnd= -1;
if ($USE_SSTFND eq 'yes'){ $use_sstfnd = 1;}
elsif ($USE_SSTFND eq 'no'){ $use_sstfnd = 0;}
else {
   die "($prog, $func) USE_SSTFND selection needs to be yes or no in configuration file ($config_file), not a valid configuration file so exiting..";
}




 # interprets choice of k parameterisation
my $k_parameterisation=0;
if ($K_HO2006 eq 'yes'){ $k_parameterisation = 1;}
elsif ($K_NIGHTINGALE2000 eq 'yes'){ $k_parameterisation = 2;}
elsif ($KT_OCEANFLUXGHG eq 'yes'){ $k_parameterisation = 3;}
elsif ($K_WANNINKHOF1992 eq 'yes'){ $k_parameterisation = 4;}
elsif ($K_MCGILLIS2001 eq 'yes'){ $k_parameterisation = 5;}
elsif ($K_HO1997 eq 'yes'){ $k_parameterisation = 6;}
elsif ($KD_OCEANFLUXGHG_BACKSCATTER eq 'yes'){ $k_parameterisation = 7;}
elsif ($KB_OCEANFLUXGHG eq 'yes'){ $k_parameterisation = 8;}
elsif ($K_GENERIC eq 'yes'){ $k_parameterisation = 9;}
elsif ($KD_OCEANFLUXGHG_WIND eq 'yes'){ $k_parameterisation = 10;}
elsif ($KT_OCEANFLUXGHG_KD_WIND eq 'yes'){ $k_parameterisation = 11;}
elsif ($K_WANNINKHOF2013 eq 'yes'){ $k_parameterisation = 12;}
elsif ($RAIN_WET_DEPOSITION eq 'yes'){ $k_parameterisation = 0;}
else {
   die "($prog, $func) No k parameterisation selected in configuration file ($config_file), not a valid configuration file so exiting..";
}

if ($debug ==1){
   print "\n\n($prog, $func) k_paramterisation: $k_parameterisation\n\n";
}

 # interprets rain controls
my $rain_sst_switch=0;
if ($BIAS_SSTSKIN_DUE_RAIN eq 'yes'){ $rain_sst_switch = 1;}
elsif ($BIAS_SSTSKIN_DUE_RAIN eq 'no') {$rain_sst_switch = 0;}
else{
   die "($prog, $func) unknown entry for bias_sstskin_due_rain selected in configuration file ($config_file), not a valid configuration file so exiting..";
}

my $rain_wet_deposition_switch=0;
if ($RAIN_WET_DEPOSITION eq 'yes'){ $rain_wet_deposition_switch = 1;}
elsif ($RAIN_WET_DEPOSITION eq 'no') {$rain_wet_deposition_switch = 0;}
else{
   die "($prog, $func) unknown entry for rain_wet_deposition selected in configuration file ($config_file), not a valid configuration file so exiting..";
}

my $k_rain_linear_ho1997_switch=0;
if ($K_RAIN_LINEAR_HO1997 eq 'yes'){ $k_rain_linear_ho1997_switch = 1;}
elsif ($K_RAIN_LINEAR_HO1997 eq 'no') {$k_rain_linear_ho1997_switch = 0;}
else{
   die "($prog, $func) unknown entry for k_rain_linear_ho1997 selected in configuration file ($config_file), not a valid configuration file so exiting..";
}

my $k_rain_nonlinear_h2012_switch=0;
if ($K_RAIN_NONLINEAR_H2012 eq 'yes'){ $k_rain_nonlinear_h2012_switch = 1;}
elsif ($K_RAIN_NONLINEAR_H2012 eq 'no') {$k_rain_nonlinear_h2012_switch = 0;}
else{
   die "($prog, $func) unknown entry for k_rain_nonlinear_h2012 selected in configuration file ($config_file), not a valid configuration file so exiting..";
}


 # displaying entries
if ($debug ==1){
   print "SRC_HOME: $SRC_HOME\n";
   print "RAPID: $RAPID\n";
   print "EQUILIBRIUM: $EQUILIBRIUM\n";
   print "BULK: $BULK\n";
   print "SST_GRADIENTS: $SST_GRADIENTS\n";   
   print "USE_SSTSKIN: $USE_SSTSKIN\n";
   print "USE_SSTFND: $USE_SSTFND\n";  
   print "SALINE_SKIN_VALUE: $SALINE_SKIN_VALUE\n";
   print "WINDU10: $WINDU10\n";
   print "SSTSKIN: $SSTSKIN\n";
   print "SSTFND: $SSTFND\n";
   print "SIGMA0: $SIGMA0\n";
   print "SIG_WV_HT: $SIG_WV_HT\n";
   print "ICE: $ICE\n";
   print "WINDU10_PROD_NAME: $WINDU10_PROD_NAME\n";
   print "SSTSKIN_PROD_NAME: $SSTSKIN_PROD_NAME\n";
   print "SSTFND_PROD_NAME: $SSTFND_PROD_NAME\n";
   print "SIGMA0_PROD_NAME: $SIGMA0_PROD_NAME\n";
   print "SIG_WV_HT_PROD_NAME: $SIG_WV_HT_PROD_NAME\n";
   print "ICE_PROD_NAME: $ICE_PROD_NAME\n";
   print "PCO2_PROD_NAME: $PCO2_PROD_NAME\n";
   print "PCO2_SST_PROD_NAME: $PCO2_SST_PROD_NAME\n";
   print "VCO2_PROD_NAME: $VCO2_PROD_NAME\n";
   print "SALINITY_PROD_NAME: $SALINITY_PROD_NAME\n";
   print "RAIN_PROD_NAME: $RAIN_PROD_NAME\n";
   print "PRESSURE: $PRESSURE\n";
   print "PCO2: $PCO2\n";
   print "SALINITY: $SALINITY\n";
   print "RAIN: $RAIN\n";
   print "PCO2_DATA_SELECTION: $PCO2_DATA_SELECTION\n";
   print "SALINITY_DATA_SELECTION: $SALINITY_DATA_SELECTION\n";
   print "RAIN_DATA_SELECTION: $RAIN_DATA_SELECTION\n";
   print "SSTGRAD: $SSTGRAD\n";
   print "BIOLOGY: $BIOLOGY\n";
   print "SUSPENDED_PARTICLES: $SUSPENDED_PARTICLES\n";
   print "ATLANTIC_OCEAN: $ATLANTIC_OCEAN\n";
   print "PACIFIC_OCEAN: $PACIFIC_OCEAN\n";
   print "SOUTHERN_OCEAN: $SOUTHERN_OCEAN\n";
   print "INDIAN_OCEAN: $INDIAN_OCEAN\n";
   print "LONGHURST: $LONGHURST\n";
   print "RANDOM_NOISE_WINDU10: $RANDOM_NOISE_WINDU10\n";
   print "RANDOM_NOISE_SSTSKIN: $RANDOM_NOISE_SSTSKIN\n";
   print "RANDOM_NOISE_SSTFND: $RANDOM_NOISE_SSTFND\n";
   print "RANDOM_NOISE_PCO2: $RANDOM_NOISE_PCO2\n";
   print "BIAS_WINDU10: $BIAS_WINDU10\n";
   print "BIAS_SSTSKIN: $BIAS_SSTSKIN\n";
   print "BIAS_SSTFND: $BIAS_SSTFND\n";
   print "BIAS_PCO2: $BIAS_PCO2\n";
   print "BIAS_WINDU10_VALUE: $BIAS_WINDU10_VALUE\n";
   print "BIAS_SSTSKIN_VALUE: $BIAS_SSTSKIN_VALUE\n";
   print "BIAS_SSTFND_VALUE: $BIAS_SSTFND_VALUE\n";
   print "BIAS_PCO2_VALUE: $BIAS_PCO2_VALUE\n";
   print "BIAS_K: $BIAS_K\n";
   print "BIAS_K_PERCENT: $BIAS_K_PERCENT\n";
   print "BIAS_K_VALUE: $BIAS_K_VALUE\n";
   print "BIAS_K_BIOLOGY_VALUE: $BIAS_K_BIOLOGY_VALUE\n";
   print "BIAS_K_WIND_VALUE: $BIAS_K_WIND_VALUE\n";
   print "K_HO2006: $K_HO2006\n";
   print "K_NIGHTINGALE2000: $K_NIGHTINGALE2000\n";
   print "KT_OCEANFLUXGHG: $KT_OCEANFLUXGHG\n";
   print "K_WANNINKHOF1992: $K_WANNINKHOF1992\n";
   print "K_WANNINKHOF2013: $K_WANNINKHOF2013\n";
   print "K_MCGILLIS2001: $K_MCGILLIS2001\n";
   print "K_HO1997: $K_HO1997\n";
   print "KD_OCEANFLUXGHG_BACKSCATTER: $KD_OCEANFLUXGHG_BACKSCATTER\n";
   print "KD_OCEANFLUXGHG_WIND: $KD_OCEANFLUXGHG_WIND\n";
   print "KB_OCEANFLUXGHG: $KB_OCEANFLUXGHG\n";
   print "KT_OCEANFLUXGHG_KD_WIND: $KT_OCEANFLUXGHG_KD_WIND\n";
   print "K_GENERIC: $K_GENERIC\n";
   print "K_GENERIC_SC: $K_GENERIC_SC\n";   
   print "K_GENERIC_A0: $K_GENERIC_A0\n";
   print "K_GENERIC_A1: $K_GENERIC_A1\n";
   print "K_GENERIC_A2: $K_GENERIC_A2\n";
   print "K_GENERIC_A3: $K_GENERIC_A3\n";      
   print "KD_WEIGHTING: $KD_WEIGHTING\n";
   print "KB_WEIGHTING: $KB_WEIGHTING\n";
   print "KB_ASYMMETRY: $KB_ASYMMETRY\n";
   print "RAIN_WET_DEPOSITION: $RAIN_WET_DEPOSITION\n";
   print "K_RAIN_LINEAR_HO1997: $K_RAIN_LINEAR_HO1997\n";
   print "K_RAIN_NONLINEAR_H2012: $K_RAIN_NONLINEAR_H2012\n";   
   print "BIAS_SSTSKIN_DUE_RAIN: $BIAS_SSTSKIN_DUE_RAIN\n";
   print "BIAS_SSTSKIN_DUE_RAIN_VALUE: $BIAS_SSTSKIN_DUE_RAIN_VALUE\n";
   print "BIAS_SSTSKIN_DUE_RAIN_WIND: $BIAS_SSTSKIN_DUE_RAIN_WIND\n";
   print "BIAS_SSTSKIN_DUE_RAIN_INTENSITY: $BIAS_SSTSKIN_DUE_RAIN_INTENSITY\n"
}

#my $processing_time = scalar(gmtime)." GMT";
#my $processing_time = strftime("%a %b %e %H:%M:%S %Y GMT", gmtime);

my ($pro_sec,$pro_min,$pro_hour,$pro_mday,$pro_mon,$pro_year,$pro_wday,$pro_yday,$pro_isdst) = gmtime(time);
$pro_year = $pro_year + 1900;
$pro_mon += 1; # count goes from 0-11, rather an 1-12 as humans read
my @pro_weekday = ("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday");
my @pro_months = ("January", "February", "March","April", "May", "June", "July", "August", "September", "October", "November", "December");
$pro_hour = sprintf("%02d",$pro_hour); 
$pro_min = sprintf("%02d",$pro_min);
$pro_sec = sprintf("%02d",$pro_sec);
my $pro_mon_formatted = sprintf("%02d",$pro_mon);
my $pro_mday_formatted = sprintf("%02d",$pro_mday);
my $processing_time = "$pro_weekday[$pro_wday]"."_"."$pro_mday$pro_months[$pro_mon-1]$pro_year"."_"."$pro_mday_formatted/$pro_mon_formatted/$pro_year"."_at_"."$pro_hour:$pro_min:$pro_sec"."GMT";
my $process_machine_time_string = "Executing on $hostname at ".$processing_time;
print "\n($prog, $func) $process_machine_time_string";
print "\n($prog, $func) Config file: $config_file";

printf STDOUT ("\n($prog, $func) Config file read successfull, all required keys found.");


 # biological data files GlobColour data
 # TO-DO: will be removed when we start using ESA Ocean Colour CCI data
#my $bio_list = "L3m_20100101-20100131__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100201-20100228__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100301-20100331__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100401-20100430__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100501-20100531__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100601-20100630__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100701-20100731__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100801-20100831__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20100901-20100930__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20101001-20101031__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20101101-20101130__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc L3m_20101201-20101231__GLOB_100_GSM-MERMODSWF_CHL1_MO_00.nc";
#my @bio_files = split /\s+/,$bio_list;

 # biological data files from ESA Ocean Colour CCI
#my $bio_list = "ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201001-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201002-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201003-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201004-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201005-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201006-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201007-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201008-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201009-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201010-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201011-fv0.95.nc ESACCI-OC-MAPPED-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_PML_OC4v6_QAA-201012-fv0.95.nc";
#my @bio_files = split /\s+/,$bio_list;


 # about to generate the climatology
my $month_list = "jan feb mar apr may jun jul aug sep oct nov dec";
my @month_names = split /\s+/,$month_list;

#for (my $i=1; $i<=scalar(@pco2_files); $i++){
 # inerpretting the random noise switches
my $random_noise_windu10_switch=0;
my $random_noise_sstskin_switch=0;
my $random_noise_sstfnd_switch=0;
my $random_noise_pco2_switch=0;
if ($RANDOM_NOISE_WINDU10 eq "yes"){$random_noise_windu10_switch=1;}
if ($RANDOM_NOISE_SSTSKIN eq "yes"){$random_noise_sstskin_switch=1;}
if ($RANDOM_NOISE_SSTFND eq "yes"){$random_noise_sstfnd_switch=1;}
if ($RANDOM_NOISE_PCO2 eq "yes"){$random_noise_pco2_switch=1;}

# interpreting the bias switches
my $bias_windu10_switch=0;
my $bias_sstskin_switch=0;
my $bias_sstfnd_switch=0;
my $bias_pco2_switch=0;
my $bias_k_switch=0;
if ($BIAS_WINDU10 eq "yes"){$bias_windu10_switch=1;}
if ($BIAS_SSTSKIN eq "yes"){$bias_sstskin_switch=1;}
if ($BIAS_SSTFND eq "yes"){$bias_sstfnd_switch=1;}
if ($BIAS_PCO2 eq "yes"){$bias_pco2_switch=1;}
if ($BIAS_K eq "yes"){$bias_k_switch=1;}

my $bias_k_percent_switch=0;
if ($BIAS_K_PERCENT eq "yes"){$bias_k_percent_switch=1;}


 # interpreting whether or not to flip data
my $pco2_data_selection = -1;
if ($PCO2_DATA_SELECTION eq "taka"){
   $pco2_data_selection = 0;
   print "\n($prog, $func) Using Takahashi climatological pco2 data (with lat grid +90 to -90 degrees).";
   }
elsif ($PCO2_DATA_SELECTION eq "socat_pco2"){
   $pco2_data_selection = 1;
   print "\n($prog, $func) Using SOCAT climatological pco2 data (with lat grid -90 to +90 degrees).";
   }
elsif ($PCO2_DATA_SELECTION eq "socat_fco2"){
   $pco2_data_selection = 2;
   print "\n($prog, $func) Using SOCAT climatological fco2 data (with lat grid -90 to +90 degrees).";
   }
  elsif ($PCO2_DATA_SELECTION eq "insitu_pco2"){
   $pco2_data_selection = 3;
   print "\n($prog, $func) Using insitu pco2 data (with lat grid -90 to +90 degrees).";
   } 
else{
  die "($prog, $func) PCO2_DATA_SELECTION is invalid in configuration file ($config_file), not a valid configuration file so exiting..";
}

if ($pco2_data_selection == 0 or $pco2_data_selection == 3){
          # Takahashi data or in-situ data
         $PCO2 = $home_dir."/".$PCO2;
}

 # limited to 2010 for climatology runs
my $year;
for (my $y=$year_start; $y<=$year_end; $y++){
   $year = sprintf("%04d",$y);
   print "\n($prog, $func) Running year: $year";
   
   for (my $i=1; $i<=12; $i++){
   #for (my $i=1; $i<=1; $i++){ # testing, running one month only
    
       # create input filenames
       # year specific datasets
      my $windu10_file = $WINDU10."/".$year."/".$year.sprintf("%02d", $i)."_OCF-WSP-???-1M-*-*.nc";
      my @file_glob = glob("$windu10_file");
      if (scalar(@file_glob) != 1){
         die "($prog, $func) More than one file or no file found in glob - arctic or global - ($windu10_file), exiting.";
      } else {
         $windu10_file = $file_glob[0];
      }
      
      #IGA Temporaary fix
      # my $ice_file = $ICE."/".$year."/".$year.sprintf("%02d", $i)."_OCF-ICE-???-1M-*-*.nc";
      # @file_glob = glob("$ice_file");
      # if (scalar(@file_glob) != 1){
      #    die "($prog, $func) More than one file or no file found in glob ($ice_file), exiting.";
      # } else {
      #    $ice_file = $file_glob[0];
      # }
      my $ice_file = $ICE."/2010/2010".sprintf("%02d", $i)."01_OCF-ICE-???-1M-*-*.nc";
      @file_glob = glob("$ice_file");
      if (scalar(@file_glob) != 1){
         die "($prog, $func) More than one file or no file found in glob ($ice_file), exiting.";
      } else {
         $ice_file = $file_glob[0];
      }
      
      my $sstskin_file = $SSTSKIN."/".$year."/".$year.sprintf("%02d", $i)."_OCF-SST-???-1M-*-*.nc";
      @file_glob = glob("$sstskin_file");
      if (scalar(@file_glob) != 1){
         die "($prog, $func) More than one file or no file found in glob ($sstskin_file), exiting.";
      } else {
         $sstskin_file = $file_glob[0];
      }
       
       # loading the biological data
      my $bio_file = $BIOLOGY."/".$year."/*".$year.sprintf("%02d", $i)."*-fv0.95.nc";
      @file_glob = glob("$bio_file");
      if (scalar(@file_glob) != 1 and $process_layers_off==0){
         die "($prog, $func) More than one file or no file found in glob ($bio_file), exiting.";
      } elsif (scalar(@file_glob) == 0 and $process_layers_off==1) {
         $bio_file = "null"
      } else {
         $bio_file = $file_glob[0];
      }
      
      my $sstfnd_file;
      if ($use_sstfnd == 1){
         $sstfnd_file = $SSTFND."/".$year."/".$year.sprintf("%02d", $i)."_OCF-SST-???-1M-*-*.nc";
         @file_glob = glob("$sstfnd_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($sstfnd_file), exiting.";
         } else {
            $sstfnd_file = $file_glob[0];
         }
      } else {
         $sstfnd_file = $sstskin_file;
      }
      
      my $sigma0_file = $SIGMA0."/".$year."/".$year.sprintf("%02d", $i)."_OCF-SI0-???-1M-*-*.nc";
      @file_glob = glob("$sigma0_file");
      if (scalar(@file_glob) != 1){
         die "($prog, $func) More than one file or no file found in glob ($sigma0_file), exiting.";
      } else {
         $sigma0_file = $file_glob[0];
      }
      
      my $sig_wv_ht_file = $SIG_WV_HT."/".$year."/".$year.sprintf("%02d", $i)."_OCF-SSH-???-1M-*-*.nc";
      @file_glob = glob("$sig_wv_ht_file");
      if (scalar(@file_glob) != 1){
         die "($prog, $func) More than one file or no file found in glob ($sig_wv_ht_file), exiting.";
      } else {
         $sig_wv_ht_file = $file_glob[0];
      }
             
      my $pressure_file = $PRESSURE."/".$year."/".$year.sprintf("%02d", $i)."_OCF-PRE-???-1M-*-*.nc";
      @file_glob = glob("$pressure_file");
      if (scalar(@file_glob) != 1){
         die "($prog, $func) More than one file or no file found in glob ($pressure_file), exiting.";
      } else {
         $pressure_file = $file_glob[0]; 
      }
      
       # salinity options
       # either Takahashi climatology salinity or SMOS
      my $salinity_file ="null";
      my $salinity_data_selection = -1;
      if ($SALINITY_DATA_SELECTION eq "taka"){
         $salinity_data_selection = 0;
    $salinity_file = $SALINITY."/*".sprintf("%03s", $month_names[$i-1])."*.nc";
         @file_glob = glob("$salinity_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($salinity_file), exiting.";
         } else {
            $salinity_file = $file_glob[0];
         }   
      } elsif ($SALINITY_DATA_SELECTION eq "smos"){
         $salinity_data_selection = 1;
         $salinity_file = $SALINITY."/".$year."/".$year.sprintf("%02d", $i)."_POA-SSS-???-1M-*-MIRAS-SMOS.nc";
         @file_glob = glob("$salinity_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($salinity_file), exiting.";
         } else {
            $salinity_file = $file_glob[0];
         }      
      } else {
         die "($prog, $func) salinity_data_selection unrecognised ($SALINITY_DATA_SELECTION), invalid configuration file so exiting.";
      }
      
      # rain data options
      # either trmm or gpcp
      my $rain_file = "null";
      my $rain_data_selection = -1;
      if ($RAIN_DATA_SELECTION eq "gpcp"){
         $rain_data_selection = 1;
         #$rain_file = $RAIN."/".sprintf("%04d",$year)."/*".sprintf("%04d%02d",$year,$i).".nc";
         $rain_file = $RAIN."/2010/gpcp_v22.2010".sprintf("%02d",$i).".nc";#IGA fixed to 2010 for a temporary fix due to lck of data 5/5/2016
         @file_glob = glob("$rain_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($rain_file), exiting.";
         } else {
            $rain_file = $file_glob[0];
         }
      } elsif ($RAIN_DATA_SELECTION eq "trmm"){
         $rain_data_selection = 0;
    $rain_file = $RAIN."/".$year."/".$year.sprintf("%02d", $i)."_OCF-RAI-???-1M-*-*.nc";
         @file_glob = glob("$rain_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($rain_file), exiting.";
         } else {
            $rain_file = $file_glob[0]; 
         }
      } else {
         die "($prog, $func) rain_data_selection unrecognised ($RAIN_DATA_SELECTION), invalid configuration file so exiting.";
      }
      
       # fixed year and/or climatology data filename creation
      my $sstgrad_file = $SSTGRAD.sprintf("%02d", $i)."-UKMO-L4LRfnd-GLOB-v01-fv02-OSTIAGRADclim.nc";
   
       # project generated kriged SOCAT data
      my $pco2_file; 
   
      if ($pco2_data_selection == 0){
         # globing for takahashi data based on month name
         $pco2_file = $PCO2."/*".sprintf("%03s", $month_names[$i-1])."*.nc";
         @file_glob = glob("$pco2_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($pco2_file), exiting.";
         } else {
            $pco2_file = $file_glob[0];
         }
      
      } elsif (($pco2_data_selection == 1) or ($pco2_data_selection == 2)){
          # SOCAT data
          # note SOCAT data are currently in a user space, so the next line is commented out
          #$PCO2 = $home_dir."/".$PCO2;
         $pco2_file = $PCO2."/*".sprintf("%02d", $i)."_OCF-CO2-???-1M-*-*.nc";
         @file_glob = glob("$pco2_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($pco2_file), exiting.";
         } else {
            $pco2_file = $file_glob[0];
         }
      } elsif ($pco2_data_selection == 3){
          # In-situ data
          # note SOCAT data are currently in a user space, so the next line is commented out
         #$PCO2 = $home_dir."/".$PCO2;
         $pco2_file = $PCO2;
         @file_glob = glob("$pco2_file");
         if (scalar(@file_glob) != 1){
            die "($prog, $func) More than one file or no file found in glob ($pco2_file), exiting.";
         } else {
            $pco2_file = $file_glob[0];
         }
            
      } else {
         die "($prog, $func) pco2_data_flip and PCO2_DATA_SELECTION is invalid in configuration file ($config_file), not a valid configuration file so exiting..";  
      }
                 
       # create output filename
      my $out_file = "OceanFluxGHG-month".sprintf("%02d", $i)."-".$month_names[$i-1]."-".$year."-v0.nc";
      
      my $final_output_dir="/null";
      if ($output_dir_override){
          # creatin output path based on year and month if output dir is set on the command line
         my $dir_month = sprintf("%02d", $i);
    my $dir_year = sprintf("%04d", $y);
    $final_output_dir = "$OUTPUT_DIR/$dir_year/$dir_month/";
          # create the output directory
    unless(-e $final_output_dir or mkpath $final_output_dir){
       die "unable to create output directory ($final_output_dir), exiting.";
    }
    
    
      } else {
           # default option is no year/month structure
      # directory assumed to exist
         $final_output_dir = $OUTPUT_DIR;
      }
       # globcolur variable name CHL1_mean CHL1_mean
       # ESA OC CCI variable names chlor_a Rrs_555 $BIOLOGY/$bio_files[$i-1] $SUSPENDED_PARTICLES/$bio_files[$i-1]
       
       # creat the command to run       
      
      my $command = $home_dir."/".$SRC_HOME."/ofluxghg-flux-calc.py"." $sstskin_file $sstfnd_file $windu10_file $pressure_file $pco2_file $salinity_file $sigma0_file $sig_wv_ht_file $rain_file $bio_file $bio_file $ice_file $ATLANTIC_OCEAN $PACIFIC_OCEAN $SOUTHERN_OCEAN $INDIAN_OCEAN $LONGHURST $sstgrad_file $final_output_dir/$out_file $SSTSKIN_PROD_NAME $SSTFND_PROD_NAME $WINDU10_PROD_NAME $SIGMA0_PROD_NAME $SIG_WV_HT_PROD_NAME $ICE_PROD_NAME $PCO2_PROD_NAME $PCO2_SST_PROD_NAME $VCO2_PROD_NAME msl_mean $SALINITY_PROD_NAME $RAIN_PROD_NAME chlor_a Rrs_555 sea-mask sea-mask gradient_fields $year $random_noise_windu10_switch $random_noise_sstskin_switch $random_noise_sstfnd_switch $random_noise_pco2_switch $bias_windu10_switch $bias_sstskin_switch $bias_sstfnd_switch $bias_pco2_switch $bias_k_switch $BIAS_WINDU10_VALUE $BIAS_SSTSKIN_VALUE $BIAS_SSTFND_VALUE $BIAS_PCO2_VALUE $BIAS_K_VALUE $bias_k_percent_switch $BIAS_K_BIOLOGY_VALUE $BIAS_K_WIND_VALUE $pco2_data_selection $salinity_data_selection $rain_data_selection $k_parameterisation $K_GENERIC_SC $K_GENERIC_A0 $K_GENERIC_A1 $K_GENERIC_A2 $K_GENERIC_A3 $KD_WEIGHTING $KB_WEIGHTING $KB_ASYMMETRY $rain_sst_switch $BIAS_SSTSKIN_DUE_RAIN_VALUE $BIAS_SSTSKIN_DUE_RAIN_INTENSITY $BIAS_SSTSKIN_DUE_RAIN_WIND $rain_wet_deposition_switch $k_rain_linear_ho1997_switch $k_rain_nonlinear_h2012_switch $flux_model $sst_gradients $use_sstskin $use_sstfnd $SALINE_SKIN_VALUE $process_layers_off $config_file $hostname $processing_time";
      print "Command dir: $command\n";
      print "END COMMAND--------\n";
      if ($debug == 1){
         print "$command \n";
      }
      
      print "\n($prog, $func) Output file: $final_output_dir/$out_file\n";
    
       # execute the command
      system("$command");
    
       # interpret the return status from the python
     if ( $? == -1 ){
        print "\n($prog, $func) Command failed with $!\n";
     } elsif ( ($? >> 8) != 0) {
        printf "\n($prog, $func) Run completed with errors (return: %d)\n", $? >> 8;
     } else {
        printf "\n($prog, $func) Run completed with no errors (return: %d)\n", $? >> 8;
     }
     
   }
}




####################################################
# subroutines 
####################################################


# Print usage information;
sub usage() {
   print <<EOF

 Util for generating the ESA OceanFlux Greenhouse Gases global climatology and time series fluxes.
 (development, use with care)
 
 usage:
    Ofluxghg-run-climatology.pl --config <onfiguration file> [-pco2_dir_override|p <input dir> -output_dir_override|o <output dir> -process_layers_off|l -year_start|s <year> -year_end|e <year>]

  year default is 2010 (ie -year_start 2010 -year_end 2010)

 AUTHOR
   Jamie Shutler <jams\@pml.ac.uk> 2012-2013
   
EOF
}
