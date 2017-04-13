library(googleComputeEngineR)
library(googleCloudStorageR)
library(googleAuthR)

gce_global_project("my-project")
gce_global_zone("europe-west1-b")
gcs_global_bucket("your-gcs-bucket")

## auth to same project we're on
googleAuthR::gar_gce_auth()

## launch the premade VM
vm <- gce_vm("slave-1")

## set SSH to use 'master' username as configured before
vm <- gce_ssh_setup(vm, username = "master", ssh_overwrite = TRUE)

## run the script on the VM that will source from GCS
runme <- "Rscript -e \"googleAuthR::gar_gce_auth();googleCloudStorageR::gcs_source('download.R', bucket = 'your-gcs-bucket')\""
out <- docker_cmd(vm, 
                  cmd = "exec", 
                  args = c("rstudio", runme), 
                  wait = TRUE)