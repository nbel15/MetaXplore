 MetaXplore

# Overview
MetaXplore, an interactive, user-friendly platform that enables the discovery and visualization of amplicon sequencing data. Currently, it provides a set of well-documented choices for downstream analysis, including alpha and beta di-versity analysis, taxonomic composition, differential abundance analysis, identification of the core microbiome within a population and biomarker analysis. These features are presented in a user-friendly format that facilitates easy customization and the generation of publication-quality graphics. MetaXplore is implemented entirely in the R language using the Shiny framework. It can be easily used locally on any system with R installed, including Windows, Mac OS, and most Linux distributions, or remotely through a web server (<a href="http://metaxplore.eu"><b>here</a>) without bioinformatic expertise. It can also be used as a framework for advanced users who can modify and expand the tool.

You can easily access MetaXplore here: http://metaxplore.eu


<p align="center">
  <img src="www/img/overview.png" height="60%" width="60%"></center> 
</p>

# Description of the content
**Folders:
 - "analysis_functions": Contains the main functions of the downstream analysis (relative abundance, Alpha diversity, ....)
 - "server_functions": Contains the server logic required for the analysis sections
 - "www": Contains the guide sections in HTML format and the example files

 **Files:
 - server.R: contains the server function definition, it is require the files saved in "server_functions" folder
 - ui.R: contains the user interface definition
 - tools_functions: contains six extra functions needed during the analysis
	 
