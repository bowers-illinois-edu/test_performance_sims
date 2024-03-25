## A file to setup clusters in different hardware environments

library(here)
## You will need to install Rmpi, and to do this on OS X you'll need to  do brew install open-mpi
## You will also need to  install snow to use Rmpi
library(parallel)
library(parallelly)
library(future)
library(future.apply)

## Get number of cores perhaps from an environment variable
if (Sys.getenv("CORES") == "" & !exists("numcores")) {
        numcores <- detectCores() ## as.numeric(system("sysctl hw.ncpu | awk '{print $2}'",intern=TRUE))
} else {
        if (!exists("numcores")) {
                numcores <- as.numeric(Sys.getenv("CORES")[[1]])
        }
}

## Different kinds of clusters
if (!exists("clusterswitch")) {
        clusterswitch <- Sys.getenv("CLUSTER")
}
stopifnot(clusterswitch %in% c("", "none", "ec2", "dkh", "taub", "keeling", "keeling-future", "stout", "keeling-socket"))

if (clusterswitch == "none") {
        ## thecluster<-rep("localhost",numcores)
        ## cl<-makeCluster(thecluster)
        cl <- NULL
}
if (clusterswitch == "") {
        ## By default use the local machine
        # thecluster <- rep("localhost", numcores)
        # cl <- makeCluster(thecluster)
        plan(multicore, workers = numcores)
}
if (clusterswitch == "ec2") {
        ## This may or may not work on Amazon EC2 instances. It once did...
        nodes <- system("grep node[0-9] /etc/hosts | cut -d ' ' -f 2", intern = TRUE)
        mastercores <- length(system("ssh master  cat /proc/cpuinfo | grep processor | cut -d ':' -f 2", intern = TRUE))
        nodecores <- length(system("ssh node001 cat /proc/cpuinfo | grep processor | cut -d ':' -f 2", intern = TRUE))
        thecluster <- c(rep("master", mastercores), rep(nodes, each = nodecores))
        cl <- makeCluster(thecluster)
}
if (clusterswitch == "dkh") {
        ## Trying to use both machines in my office
        nodes <- c(rep("localhost", numcores), rep("jwbowers.pol.illinois.edu", 8))
        thecluster <- nodes
        cl <- makeCluster(thecluster)
}
if (clusterswitch == "taub" & Sys.getenv("PBS_NODEFILE") != "") {
        ## This is for the campus cluster not keeling
        thecluster <- system("cat $PBS_NODEFILE|wc -l", intern = TRUE)
        library(Rmpi)
        cl <- makeCluster(thecluster[1], type = "MPI")
}
if (clusterswitch == "keeling") {
        library(Rmpi)
        cl <- makeCluster(mpi.universe.size(), type = "MPI")
        ## 	plan(cluster,workers=cl)
}
if (clusterswitch == "keeling-future") {
	workers <- parallelly::availableWorkers()
	# nworkers <- length(workers)
        cl <- parallelly::makeClusterPSOCK(workers,verbose=TRUE)
	## Errors about non-exportable objects from data.table are false positives
	## https://cran.r-project.org/web/packages/future/vignettes/future-4-non-exportable-objects.html
	plan(cluster, workers = cl,persistent=TRUE)
}
if (clusterswitch == "keeling-socket") {
        workers <- parallelly::availableWorkers()
        cl <- parallelly::makeClusterPSOCK(workers, verbose = TRUE,autoStop = TRUE)
}
