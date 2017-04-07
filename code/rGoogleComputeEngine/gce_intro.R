

library(googleComputeEngineR)
library(googleCloudStorageR)
library(googleAuthR)
setwd(file.path("/Users",
                "taylerblake",
                "Documents",
                "Dissertation",
                "code"))

Sys.setenv(GCE_SSH_USER="Tayler Blake")
Sys.setenv(GCE_AUTH_FILE = file.path("/Users",
                                     "taylerblake",
                                     "Documents",
                                     "Dissertation",
                                     "code",
                                     "rGoogleComputeEngine",
                                     "Dissertation Computation-ac8e343f99cc.json"))

## set options for authentication
options(googleAuthR.scopes.selected = c("https://www.googleapis.com/auth/cloud-platform",
                                        "https://www.googleapis.com/auth/compute"))

gce_global_project(project = "linen-creek-99303")
gce_global_zone("us-east1-b")
## authenticate
## using service account, ensure service account email added to GA account, BigQuery user permissions set, etc.
googleAuthR::gar_auth_service(file.path("/Users",
                                        "taylerblake",
                                        "Documents",
                                        "Dissertation",
                                        "code",
                                        "rGoogleComputeEngine",
                                        "Dissertation Computation-ac8e343f99cc.json"))

r_vm <- gce_vm_template(template = "r-base",
             name = "r-vm",
             predefined_type = "n1-highmem-2")


docker_build(r_vm, 
             dockerfolder = file.path("/Users",
                                    "taylerblake",
                                    "Documents",
                                    "Dissertation",
                                    "code",
                                    "rGoogleComputeEngine",
                                    "dockerFolder"), 
             new_image = "my-r-image")

gce_push_registry(r_vm, "my-r-image", image_name = "my-r-image")
vm3 <- gce_vm(template = "rstudio",
              name = "my-rstudio-3",
              username = "mark", password = "mark1234",
              predefined_type = "n1-highmem-2",
              dynamic_image = gce_tag_container("my-custom-image"))

gcs_create_bucket("dissertation-computation-bucket",
                  projectId=gce_get_project()$name, 
                  location = "US",
                  storageClass = c("STANDARD"),
                  predefinedAcl = c("private"),
                  predefinedDefaultObjectAcl = c("private"),
                  projection = c("noAcl","full"),
                  versioning = TRUE,
                  lifecycle = NULL)

gcs_upload(,
           "dissertation-computation-bucket")

gce_ssh_browser(r_vm)





## upload to Google Cloud Storage
gcs_upload(gadata, name = paste0("gadata_",Sys.Date(),".csv"))
gce_ssh_browser("my-server")




