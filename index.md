# R Interface

## Customizing Startup

You can customize the R environment through a site initialization file or a directory initialization file. R will always source the Rprofile.site file first. On Windows, the file is in the C:\Program Files\R\R-n.n.n\etc directory. You can also place a .Rprofile file in any directory that you are going to run R from or in the user home directory.

At startup, R will source the Rprofile.site file. It will then look for a .Rprofile file to source in the current working directory. If it doesn't find it, it will look for one in the user's home directory. There are two special functions you can place in these files. .First( ) will be run at the start of the R session and .Last( ) will be run at the end of the session.

```
# Sample Rprofile.site file 

# Things you might want to change
# options(papersize="a4") 
# options(editor="notepad") 
# options(pager="internal")

# R interactive prompt 
# options(prompt="> ")
# options(continue="+ ") 

# to prefer Compiled HTML 
help options(chmhelp=TRUE) 
# to prefer HTML help 
# options(htmlhelp=TRUE) 

# General options 
options(tab.width = 2) 
options(width = 130)
options(graphics.record=TRUE) 

.First <- function(){
 library(Hmisc)
 library(R2HTML)
 cat("\nWelcome at", date(), "\n") 
}

.Last <- function(){ 
 cat("\nGoodbye at ", date(), "\n")
}
```

## Getting Help
Once R is installed, there is a comprehensive built-in help system. At the program's command prompt you can use any of the following:

```
help.start()   # general help
help(foo)      # help about function foo
?foo           # same thing 
apropos("foo") # list all functions containing string foo
example(foo)   # show an example of function foo

# search for foo in help manuals and archived mailing lists
RSiteSearch("foo")

# get vignettes on using installed packages
vignette()      # show available vingettes
vignette("foo") # show specific vignette
```

### Sample Datasets
R comes with a number of sample datasets that you can experiment with. Type data( ) to see the available datasets. The results will depend on which packages you have loaded. Type help(datasetname) for details on a sample dataset.

## The Workspace

The workspace is your current R working environment and includes any user-defined objects (vectors, matrices, data frames, lists, functions). At the end of an R session, the user can save an image of the current workspace that is automatically reloaded the next time R is started. Commands are entered interactively at the R user prompt. Up and down arrow keys scroll through your command history. 

You will probably want to keep different projects in different physical directories. Here are some standard commands for managing your workspace.

```
getwd() # print the current working directory - cwd 
ls()    # list the objects in the current workspace

setwd(mydirectory)      # change to mydirectory
setwd("c:/docs/mydir")  # note / instead of \ in windows 
setwd("/usr/rob/mydir") # on linux

# view and set options for the session
help(options) # learn about available options
options() # view current option settings
options(digits=3) # number of digits to print on output

# work with your previous commands
history() # display last 25 commands
history(max.show=Inf) # display all previous commands

# save your command history 
savehistory(file="myfile") # default is ".Rhistory" 

# recall your command history 
loadhistory(file="myfile") # default is ".Rhistory"

# save the workspace to the file .RData in the cwd 
save.image()

# save specific objects to a file
# if you don't specify the path, the cwd is assumed 
save(object list,file="myfile.RData")

# load a workspace into the current session
# if you don't specify the path, the cwd is assumed 
load("myfile.RData")

q() # quit R. You will be prompted to save the workspace.

```

## Packages
Packages are collections of R functions, data, and compiled code in a well-defined format. The directory where packages are stored is called the library. R comes with a standard set of packages. Others are available for download and installation. Once installed, they have to be loaded into the session to be used.

```
.libPaths() # get library location 
library()   # see all packages installed 
search()    # see packages currently loaded
```

### Adding Packages
You can expand the types of analyses you do be adding other packages. A complete list of contributed packages is available from CRAN.

Follow these steps:

Download and install a package (you only need to do this once).
To use the package, invoke the library(package) command to load it into the current session. (You need to do this once in each session, unless you customize your environment to automatically load it each time.)
On MS Windows:

Choose Install Packages from the Packages menu.
Select a CRAN Mirror. (e.g. Norway)
Select a package. (e.g. boot)
Then use the library(package) function to load it for use. (e.g. library(boot))
On Linux:

Download the package of interest as a compressed file.
At the command prompt, install it using 
R CMD INSTALL [options] [l-lib] pkgs
Use the library(package) function within R to load it for use in the session.

### Creating Your Own Packages
To create your own packages look at Writing R Extensions (the definitive guide), Leisch's Creating R Packages: A Tutorial, and Rossi's Making R packages Under Windows: A Tutorial.

## Input/Output

By default, launching R starts an interactive session with input from the keyboard and output to the screen. However, you can have input come from a script file (a file containing R commands) and direct output to a variety of destinations.

### Input
The source( ) function runs a script in the current session. If the filename does not include a path, the file is taken from the current working directory.

```
# input a script
source("myfile")
```

### Output
The sink( ) function defines the direction of the output.

```
# direct output to a file 
sink("myfile", append=FALSE, split=FALSE)

# return output to the terminal 
sink()
```

The append option controls whether output overwrites or adds to a file. The split option determines if output is also sent to the screen as well as the output file.

Here are some examples of the sink() function.

```
# output directed to output.txt in c:\projects directory.
# output overwrites existing file. no output to terminal. 
sink("c:/projects/output.txt")

# output directed to myfile.txt in cwd. output is appended
# to existing file. output also send to terminal. 
sink("myfile.txt", append=TRUE, split=TRUE)
```

When redirecting output, use the cat( ) function to annotate the output.

### Graphs

sink( ) will not redirect graphic output. To redirect graphic output use one of the following functions. Use dev.off( ) to return output to the terminal.

```
Function	Output to
pdf("mygraph.pdf")	pdf file
win.metafile("mygraph.wmf")	windows metafile
png("mygraph.png")	png file
jpeg("mygraph.jpg")	jpeg file
bmp("mygraph.bmp")	bmp file
postscript("mygraph.ps")	postscript file
```

Use a full path in the file name to save the graph outside of the current working directory.

```
# example - output graph to jpeg file 
jpeg("c:/mygraphs/myplot.jpg")
plot(x)
dev.off()
```

## Graphic User Interfaces

R is a command line driven program. The user enters commands at the prompt ( > by default ) and each command is executed one at a time.

There have been a number of attempts to create a more graphical interface, ranging from code editors that interact with R, to full-blown GUIs that present the user with menus and dialog boxes.

RStudio is my favorite example of a code editor that interfaces with R for Windows, MacOS, and Linux platforms.

Perhaps the most stable, full-blown GUI is R Commander, which can also run under Windows, Linux, and MacOS (see the documentation for technical requirements).
https://socialsciences.mcmaster.ca/jfox/Misc/Rcmdr/

Both of these programs can make R a lot easier to use.

## Publication Quality Output

Compared with SAS and SPSS, R's ability to output results for publication quality reports is somewhat rudimentary (although this is evolving).

The R2HTML package lets you output text, tables, and graphs in HTML format. Here is a sample session, followed by an explanation.

```
# Sample Session 
library(R2HTML)
HTMLStart(outdir="c:/mydir", file="myreport",
   extension="html", echo=FALSE, HTMLframe=TRUE)
HTML.title("My Report", HR=1)

HTML.title("Description of my data", HR=3)
summary(mydata) 

HTMLhr()

HTML.title("X Y Scatter Plot", HR=2)
plot(mydata$y~mydata$x)
HTMLplot() 

HTMLStop()
```

Once you invoke HTMLStart( ), the prompt will change to HTML> until you end with HTMLStop().

The echo=TRUE option copies commands to the same file as the output.

HTMLframe=TRUE creates framed output, with commands in the left frame, linked to output in the right frame. By default, a CSS file named R2HTML.css controlling page look and feel is output to the same directory. Optionally, you can include a CSSFile= option to use your own formatting file.

Use HTML.title() to annotate the output. The HR option refers to HTML title types (H1, H2, H3, etc.) . The default is HR=2.

HTMLhr() creates a horizontal rule.

Since several interactive commands may be necessary to create a finished graph, invoke the HTMLplot() function when each graph is ready to output.

The RNews article The R2HTML Package has more complex examples using titles, annotations, header and footer files, and cascading style sheets.

### Other Options
The R Markdown Package from R Studio supports dozens of static and dynamic output formats including HTML, PDF, MS Word, scientific articles, websites, and more. (To practice R Markdown, try this tutorial taught by Garrett Grolemund, Data Scientist for R Studio.)

Sweave allows you to imbed R code in LaTeX, producing attractive reports if you know that markup language.

The odfWeave package has functions that allow you to imbedd R output in Open Document Format (ODF) files. These are the types of files created by OpenOffice software.

The SWordInstaller package allows you to add R output to Microsoft Word documents.

The R2PPT provides wrappers for adding R output to Microsoft PowerPoint presentations.

## Batch Processing

You can run R non-interactively with input from infile and send output (stdout/stderr) to another file. Here are examples.

```
# on Linux 
R CMD BATCH [options] my_script.R [outfile]

# on Microsoft Windows (adjust the path to R.exe as needed) 
"C:\Program Files\R\R-2.13.1\bin\R.exe" CMD BATCH 
   --vanilla --slave "c:\my projects\my_script.R"
```

Be sure to look at the section on I/O for help writing R scripts.

See an Introduction to R (Appendix B) for information on the command line options.
https://cran.r-project.org/doc/manuals/R-intro.pdf

## Reusing Results

In SAS, you can save the results of statistical analyses using the Output Delivery System (ODS). While ODS is a vast improvement over PROC PRINTO, it's sophistication can make some features very hard to learn (just try mastering PROC TEMPLATE). In SPSS you can do the same thing with the Output Management System (OMS). Again, not one of the easiest topics to learn.

One of the most useful design features of R is that the output of analyses can easily be saved and used as input to additional analyses.

```
# Example 1 
lm(mpg~wt, data=mtcars)
```

This will run a simple linear regression of miles per gallon on car weight using the data frame mtcars. Results are sent to the screen. Nothing is saved.

```
# Example 2 
fit <- lm(mpg~wt, data=mtcars)
```

This time, the same regression is performed but the results are saved under the name fit. No output is sent to the screen. However, you now can manipulate the results.

```
# Example 2 (continued...)

str(fit) # view the contents/structure of "fit"
```

The assignment has actually created a list called "fit" that contains a wide range of information (including the predicted values, residuals, coefficients, and more.

```
# Example 2 (continued again)
# plot residuals by fitted values
plot(fit$residuals, fit$fitted.values)
```

To see what a function returns, look at the value section of the online help for that function. Here we would look at help(lm).

The results can also be used by a wide range of other functions.

```
# Example 2 (one last time, I promise)

# produce diagnostic plots
plot(fit) 

# predict mpg from wt in a new set of data 
predict(fit, mynewdata)

# get and save influence statistics 
cook <- cooks.distance(fit)
```
