### Module description

#### Introduction
In order to create a workable module system, we need to set some ground rules 
first. Here, I will list all required features.

The most important part of a successful module system is that modules are easily
modified and swapped. In other words, it's crucial to have consistent input / 
output options.  

#### Module - parsing
For the parsing module, the following rules have been determined:

* Input is a fastq(.gz) file. 
* Input options can alter output, but file-naming should be as consistent as
	possible.
* Output is a directory for all output files, with the following naming:
	* gDNA									-		[basename]_gDNA.fastq.gz
	* statistics						-		[basename]_statistics.txt
	* [Additional files]		-		[basename]_(VARIABLE).xxx
* Output basename can be changed by input.
* Output can be stdout, but the accompanying files (statistics) will still
	be written.
* The statistics file should contain:
	 * Number of input reads.
	 * Number of filtered reads.
	 * Number of hits for each element.