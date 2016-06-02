/* This program contains standard definitions of values used in programs */

# include "bioplib/SysDefs.h"

/* Maximum filename length */

# define MAX_FILENAME_LENGTH 2000

/* abnum */

# define PROGRAM_CONFIG_FILE "/home/bsm2/abhi/DATA/Program_Config.txt"
# define MAXBUFF 	240
# define MAXSEQ  	32
# define MAXSECTIONS 	14
# define MAXRESINSECT 	40
# define LABELLENGTH	8
# define LIGHT_FV_LEN	129
# define HEAVY_FV_LEN    149 /* These should be updated if the files */
# define MAXLABEL        140 /* Must be > MAX(LIGHT_FV_LEN,HEAVY_FV_LEN) */
# define CONSENSUS_FILE  "AbConsensus.pir"
# define NUMBERING_FILE  "KabatNumbering.dat"
# define CONSENSUS_VAR   "KABATALIGN"
# define MUTATION_MATRIX "consensus.mat"
# define MODE_NORMAL     0
# define MODE_REVERSED   1
# define TRUE 1
# define FALSE 0
# define TEMPDIR "/tmp"

/* chaintype */

# define KAPPA 1
# define LAMBDA 2
# define HEAVY 3
# define ANTIGEN 4
# define SSEARCH33PATH "/acrm/usr/local/bin"
# define SUBIM_PATH "/home/bsm/martin/bin"
# define HUMAN_HEAVY_DB_PATH "/acrm/www/html/abhi/Z_SCORE_ANALYSIS/HUMAN_DATABASES/human_heavy.db.out"
# define HUMAN_LAMBDA_DB_PATH "/acrm/www/html/abhi/Z_SCORE_ANALYSIS/HUMAN_DATABASES/human_lambda.db.out"
# define HUMAN_KAPPA_DB_PATH "/acrm/www/html/abhi/Z_SCORE_ANALYSIS/HUMAN_DATABASES/human_kappa.db.out"
# define DATABASE_MAP_FILE "/acrm/www/cgi-bin/abhi/Z_SCORE_ANALYSIS/database_map_file"
# define MINIMUM_LENGTH_OF_ALIGNMENT 94
# define kappaThresholdZScore -3.873
# define lambdaThresholdZScore -4.497
# define heavyThresholdZScore -3.063

/* ga */

# define MAX_POPULATION 30000
# define MAX_NUMBER_OF_INTERFACE_POSITIONS 200
# define NUMBER_OF_FOLDS 5
# define MAX_NUMBER_OF_FOLDS 20
# define MAX_NUMBER_OF_PDB 1000
# define TEMPORARY_FILES_DIRECTORY "/acrm/home/abhi/tmp/SNNS/"
# define NUMBER_OF_HIDDEN_NODES_DEFAULT 10
# define MUTATION_RATE_DEFAULT 0.001
# define POSITIONS_TO_SWAP_DEFAULT 10

/* intermediate */

# define NUMBER_OF_RESIDUE_PROPERTIES_IN_MATRIX 4

/* snns */

# define SNNS_PATH_DEFAULT "/acrm/usr/local/apps/snns/SNNSv4.2/bin"
# define CORRELATION_PROGRAM_PATH "/home/bsm/martin/bin/correlation"
# define NUMBER_OF_SNNS_TRAINING_PARAMETERS 9
# define NUMBER_OF_RESIDUE_PARAMETERS 4
# define MAX_NUMBER_OF_INTERFACE_POSITIONS 200
# define NUMBER_OF_BINS_DEFAULT 20

/* subim */

# define MAX_NB_SEQ 26
# define MAX_SEQ_GEN 21
