

library(googleComputeEngineR)
library(googleCloudStorageR)



Sys.setenv(GCE_DEFAULT_PROJECT_ID="linen-creek-99303")
Sys.setenv(GCE_GLOBAL_PROJECT_ID="linen-creek-99303")
Sys.setenv(GCE_DEFAULT_ZONE="us-east1-b")
Sys.setenv(GCE_SSH_USER="Tayler Blake")
Sys.setenv(GCE_AUTH_FILE =file.path("C:/Users/tblake/downloads","Dissertation Computation-fdb6d93d9d86.json"))

## set options for authentication
options(googleAuthR.scopes.selected = c("https://www.googleapis.com/auth/cloud-platform",
                                        "https://www.googleapis.com/auth/compute"))


gce_auth(file.path("C:/Users/tblake/downloads","Dissertation Computation-fdb6d93d9d86.json"))
gce_global_project(project = "linen-creek-99303")
gce_global_zone("us-east1-b")
## authenticate
## using service account, ensure service account email added to GA account, BigQuery user permissions set, etc.
googleAuthR::gar_auth_service(file.path("C:/Users/tblake/downloads","Dissertation Computation-fdb6d93d9d86.json"))



r_vm <- gce_vm(template = "r-base",
             name = "r-vm",
             predefined_type = "n1-highmem-2")

gce_ssh_browser(r_vm)





gcs_create_bucket("linen-creek-bucket", projectId=gce_get_project()$id, 
                  location = "US",
                  storageClass = c("STANDARD"),
                  predefinedAcl = c("private"),
                  predefinedDefaultObjectAcl = c("private"),
                  projection = c("noAcl","full"), versioning = FALSE, lifecycle = NULL)


## get Google Analytics data
gadata <- google_analytics_4(123456, 
                             date_range = c(Sys.Date() - 2, Sys.Date() - 1),
                             metrics = "sessions",
                             dimensions = "medium",
                             anti_sample = TRUE)

## upload to Google BigQuery
bqr_upload_data(projectId = "myprojectId", 
                datasetId = "mydataset",
                tableId = paste0("gadata_",format(Sys.Date(),"%Y%m%d")),
                upload_data = gadata,
                create = TRUE)

## upload to Google Cloud Storage
gcs_upload(gadata, name = paste0("gadata_",Sys.Date(),".csv"))

gce_ssh_browser("my-server")
get_dockerfile("hadleyverse-crontab")




