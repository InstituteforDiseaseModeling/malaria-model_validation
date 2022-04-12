library(stringr)
# library(DBI)
filepath = 'C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder/data/Nigeria/Garki/garki.sql'
# 
# conn = dbConnect()
# 
# df <- dbGetQuery(conn, statement = read_file(filepath))
# 
# dbDisconnect(con)
# close(con)
# 
# 
# 
# getSQL <- function(filepath)
# {
#   con = file(filepath, "r")
#   sql.string <- ""
# 
#   while ( TRUE )
#   {
#     line <- readLines(con, n = 1)
#     if ( length(line) == 0 ) { break }
# 
#     line <- gsub("\\t", " ", line)
# 
#     if(grepl("--",line) == TRUE)
#     {
#       line <- paste(sub("--","/*",line),"*/")
#     }
# 
#     sql.string <- paste(sql.string, line)
#   }
# 
#   close(con)
#   return(sql.string)
# }
# 
# 
# garki = getSQL(filepath)
# 
# 
# 
# install.packages(c("dbplyr", "RSQLite"))
# library(dplyr)
# library(dbplyr)
# 
# test_gark <- DBI::dbConnect(RSQLite::SQLite(), filepath)



filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/garki_sql_parademo_subset.txt"
garki_para = file(filepath_ref, 'r')


get_df_from_sql <- function(filepath_ref)
{
  print('Converting formats the long way... this may take a while.')
  con = file(filepath_ref, "r")
  on_header = FALSE
  on_main = FALSE
  column_names = c()
  count_line = 0

  while ( TRUE ){
    count_line = count_line+1
    if((count_line %% 10000) == 0) print(paste0('On line ', count_line, ' out of ', file.info(filepath_ref)$size))
    
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) { break }

    line <- gsub("\\t", " ", line)

    if(grepl("--",line)){
      line <- paste(sub("--","/*",line),"*/")
    }
    
    # get the column headers
    if(on_header){
      if (grepl('TYPE=', line)){
        on_header=FALSE
        on_main = TRUE
        
        # create data frame
        df = data.frame(matrix(ncol=length(column_names),nrow=0, dimnames=list(NULL, column_names)))
        df_0 = data.frame(matrix(NA, ncol=length(column_names),nrow=1, dimnames=list(NULL, column_names)))
        # break
      } else{
        column_names = c(column_names, str_extract(line, '\\w+'))
      }
    }
    if(grepl("CREATE TABLE", line)){
      on_header = TRUE
    }
    
    # get the main entries
    if (on_main){
      if (grepl('INSERT INTO', line)){
        values = strsplit(sub("\\).*", "", sub(".*\\(", "", line)), ',')[[1]]
        df_0[1,] =values
        df = rbind(df, df_0)
      }
    }
  }
  close(con)
  return(df)
}


garki = getSQL(filepath)