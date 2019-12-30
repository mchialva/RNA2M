library(stringr)
# This script generates bash files containing commands to process raw short Illumina reads (50 bp, SE)
# to perform taxonomical assignments as described in Chialva et al. (2020).
#
# Dependencies
# the script assumes that the following dependencies are available:
# Trimmomatic, STAR, SortMeRNA, Bowtie2, taxoner64
#
# set options, reference databases and parameters
# set number of threads to use
threads<-"50"
# path to raw sequences (adapter-free reads)
data_folder<-"/home/matteo/Scrivania/Work on Server/Mycoplant_MT/"
# path to trimmomatic .jar file
trimmomatic_bin<-"/home/user/trimmomatic-0.35.jar"
# path to Host genome fasta (must be a *.fa file) and annotation (must be a .gff3 or GTF file)
host_genome_path<-"/home/user/host_genome"
# path to host and human transcriptomes .fasta to remove residual host contaminations
# e.g. create a comprehensive S. lycopersicum and H. sapiens data-base using dbcreator command in taxoner (Pongor et al., 2014).
# must be a single bowtie2 index
bowtie2_hosts<-"/home/user/Host_and_contaminants_DB_index"
# path containing host ribosmal sequences and default SortMeRNA rRNA data-bases (.fasta files)
# Custom host rRNA can be retreived at NCBI using the following query:
# "Solanum lycopersicum"[Organism] AND rRNA NOT variant NOT PREDICTED NOT complete genome NOT mrna
sortmerna_rRNA<-"/home/user/sortmerna/reference"
# path to taxoner indexed database (see Create_taxoner_NCBI_index.sh to generate Bowtie2 index files)
taxoner_DB<-"/home/user/taxoner_DB/"
# path to NCBI taxonomy (nodes.dmp file)
tax_path<-"/home/user/taxonomy/nodes.dmp"
# taxoner64 LCA neighbor score
neighbor_score<-0.97

### list input files (adapter-free reads from Chialva et al. 2018)
fastq_files<-list.files(data_folder, pattern = ".fastq.gz")

######### TRIMMOMATIC #########
### TRIM READS at length >45 and quality >28 using TRIMMOMATIC
dir.create(file.path(data_folder, "1.Trimmomatic", sep=""))
TRIMMO_outputs<-str_replace(fastq_files, ".fastq.gz", "_TRIMMOMATIC.fastq.gz")
# create commands list
TRIMMO_cmd<-paste("\njava -jar", trimmomatic_bin, "SE", "-threads", threads, CTDPT.files,
                  TRIMMO_outputs ,
                  "SLIDINGWINDOW:10:28 MINLEN:45\n")
# export bash file
cat(file=file.path(data_folder, "1.TRIMMOMATIC.sh"), TRIMMO_cmd, append=F)

######### STAR - Tomato host genome #########
#Build STAR Index files
STAR_index_cmd<-paste("\nSTAR --runMode genomeGenerate --runThreadN", threads,
                      "--genomeDir", file.path(host_genome_path, "STAR/"),
                      "--genomeFastaFiles", list.files(host_genome_path, pattern = ".fasta|.fa", full.names = T),
                      "--sjdbOverhang 100",
                      "--sjdbGTFfile", list.files(host_genome_path, pattern = ".gtf|.gff3", full.names = T),
                      "--sjdbGTFtagExonParentTranscript Parent\n")

# align cleaned reads on host genome using STAR
# adjust parameters according to intron size ranges
dir.create(paste(data_folder, "2.STAR", sep=""))
STAR_outputs<-paste(str_replace(str_replace(TRIMMO_outputs, ".fastq.gz", ""), "2.Trimmomatic", "3.STAR"), "_noHost_", sep="")
STAR_cmd<-paste("\nSTAR --runMode alignReads --alignIntronMin 40 --alignIntronMax 23000",
                "--outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1000",
                "--alignEndsType EndToEnd --runThreadN", threads,
                "--outFileNamePrefix", STAR_outputs,
                "--genomeDir", file.path(host_genome_path, "STAR/"),
                "--readFilesIn", TRIMMO_outputs, "--readFilesCommand zcat\n")

# export bash file
cat(file=file.path(data_folder, "2.STAR.sh"), STAR_index_cmd, STAR_cmd, append=F)

######### Remove residual host reads - Tomato SL.2.5.0 ########
#Index database with Bowtie2-build
#system(bowtie2-build Host_SLY_HSA.fasta Host_SLY_HSA --large-index --threads 50)
#bowtie2 -x /home/mchialva/DB/taxoner/bowtie2_update/HostDB/Host_SLY_HSA.fasta /home/mchialva/Data/Mycoplant_MT/3.STAR_SL2.50/Sample_BAL-R1_R1_CTDPT_TRIMMOMATIC_SL2.50_Unmapped.out.mate1 --un /home/mchialva/Data/Mycoplant_MT/unmapped.fq -D 25 -R 4 -N 1 -L 15 -i S,1,0.50 -p 50

# Align unmapped reads from previous step
dir.create(paste(data_folder, "3.BOWTIE2_hosts", sep=""))

# set input and outputs
STAR_unmapped<-paste(STAR_outputs, "Unmapped.out.mate1", sep="")   
bowtie2_outputs<-paste(str_replace(str_replace(STAR_unmapped, "Unmapped.out.mate1", ""), "3.STAR", "4.BOWTIE2_hosts"), "cleaned.fastq", sep="")

# rename unmapped reads as .fastq
change_name<-paste("\nmv", STAR_unmapped, paste(STAR_unmapped, ".fastq", sep=""),"\n")
# align on host and contaminants sequences to clean reads
BOWTIE2_cmd<- paste("\nbowtie2 -x", bowtie2_hosts, paste(STAR_unmapped, ".fastq", sep=""), "--un", bowtie2_outputs,
                    "-D 25 -R 4 -N 1 -L 15 -i S,1,0.50 -p", threads, ">", 
                    str_replace(bowtie2_outputs, "cleaned.fastq","host_aligned.sam\n"))
# export bash file
cat(file=file.path(data_folder, "3.Bowtie2.sh"), change_name, BOWTIE2_cmd, append=F)

######### SortMeRNA - DB Indexing #########
#### Creating rRNA indexes ###
LRNA<- " -L 14" #-- index seed length
# list .fasta and create indexes
rRNA_DB<- list.files(sortmerna_rRNA, pattern=".fasta$", full.names = T)
# create fasta-index list
indexdb_rna=NULL
for (i in 1:length(rRNA_DB)){
  seq<-rRNA_DB[i]
  db<-paste(seq, ",", paste(str_replace(substr(seq,1, str_length(seq)-6), "reference", "indexes"), "_L", substr(LRNA,5,8), sep=""), sep="")
  indexdb_rna<-paste(indexdb_rna, db, sep=":")
}
# index rRNA .fasta files
indexdb_rna_cmd<- paste("\nindexdb_rna --ref ", substr(indexdb_rna,2, str_length(indexdb_rna)), paste(LRNA, "-m 50000"), sep="") 
# export bash file
cat(file=paste(data_folder, "4.Sortmerna_indexes.sh", sep=""), indexdb_rna_cmd, append=F)

######### SortMeRNA - Ribosomal reads filtering #########
dir.create(paste(data_folder, "4.Sortmerna", sep="")) 
sortmerna_outputs<- paste(str_replace(str_replace(bowtie2_outputs, "3.BOWTIE2_hosts", "4.Sortmerna"),".fastq", "_rRNA_free"), sep="")

sortmerna_cmd<-paste("\n","sortmerna --ref ", substr(indexdb_rna,2, str_length(indexdb_rna)),
                      " --reads ", bowtie2_outputs, " --aligned ",
                      paste(str_replace(sortmerna_outputs, "_free", " "), " --other ", sortmerna_outputs,
                      " --fastx", " --blast 1", " --log", " -a ", threads," -v",
                      " --passes 14,7,1", " -e 1\n", sep=""), sep="")
# export bash file
cat(file=paste(data_folder, "4.Sortmerna.sh", sep=""), sortmerna_cmd, append=F)

# list rRNA-free reads
Sortmerna_unmapped<- paste(str_replace(str_replace(bowtie2_outputs, "3.BOWTIE2_hosts", "4.Sortmerna"),".fastq", "_rRNA_free.fastq"), sep="")

######## TAXONER #######
#create bowtie2 configuration file for taxoner64
cat(file=paste(data_folder, "extra_commands.txt", sep=""), "-D 25 -R 4 -N 1 -L 15 -i S,1,0.50")

#### run taxoner64 ###
dir.create(paste(data_folder, "5.Taxoner64", sep="")) 

taxoner64_cmd<-paste("\ntaxoner64 --dbPath ", taxoner_DB,
                     " --largeGenome --taxpath ", tax_path,
                     " --seq ", Sortmerna_unmapped,
                     " --output ", str_replace(str_replace(Sortmerna_unmapped, "4.Sortmerna", "5.Taxoner64"), ".fastq", "_taxoner64"),
                     " --bwt2-params ", paste(data_folder, "extra_commands.txt", sep=""), " --neighbor-score ", neighbor_score,
                     " --megan --threads ", threads, "\n", sep="")
# export bash file
cat(file=paste(data_folder, "5.Taxoner64.sh", sep=""), taxoner64_cmd, append=F)