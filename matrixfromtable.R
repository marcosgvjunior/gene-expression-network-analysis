#Institute of Physics - Federal University of Rio de Janeiro
#Interdisciplinary Academic Master's Degree in Applied Physics
#Student: Marcos Vieira
#October, 2019

#Instituto Oswaldo Cruz (PGBCS/IOC/Fiocruz)
#Ph.D. program on Computational and Systems Biology
#Student: Marcos Vieira
#September, 2022

library( dplyr )
library( readxl )
library( writexl )

matrixfromtable <- function( filename, actionNumber, outputName ){

  options( warn = -1 )

  storeoutputs          <- list()

  # Import and store data
  network               <- data.frame( read_excel( filename ) )
  nodes                 <- data.frame( network %>% select( nodeout, nodein, action ) )
  names( nodes )        <- c( "nodeout", "nodein", "action" )

  nnames                <- dplyr::union( nodes[[1]], nodes[[2]] )
  nnamessort            <- sort( nnames )

  message( "\nnodes list:" )
  message( nnamessort, "\n" )

  totalnodes            <- length( nnames )

  message( "Total number of nodes:" )
  message( totalnodes, "\n" )

  addNodes                <- as.data.frame( cbind( nnamessort, nnamessort ) )
  names( addNodes )       <- c( "nodeout", "nodein")
  addNodes                <- addNodes %>% mutate( action = as.integer( 0 ) )

  nodesaction             <-  bind_rows( nodes, addNodes )

  partialmatrix           <- as.data.frame( table( nodesaction ) ) %>% filter( action == actionNumber ) %>% mutate( action = 1 ) %>% select( nodeout, nodein, Freq )

  tabcompactvation        <- table( partialmatrix, deparse.level = 2 )

  dirmatrix               <- matrix( tabcompactvation, ncol=2*totalnodes )[, (totalnodes + 1):(2*totalnodes)]
  # adjmatrix               <- ( (dirmatrix + t(dirmatrix)) > 0 )*1

  datadirmatrix            <- data.frame( dirmatrix, row.names = nnamessort )
  names( datadirmatrix )   <- nnamessort

  # dataadjmatrix            <- data.frame( adjmatrix, row.names = nnamessort )
  # names( dataadjmatrix )   <- nnamessort

  storeoutputs$dirmatrix   <- datadirmatrix
  # storeoutputs$adjmatrix   <- dataadjmatrix
  storeoutputs$nnamessort  <- nnamessort

  write_xlsx( datadirmatrix, outputName )

  return(storeoutputs)
}

# outstore   <- matrixfromtable( "table.xlsx", 1, "act.xlsx" )
# outstore   <- matrixfromtable( "table.xlsx", 2, "sup.xlsx" )

# outstore$dirmatrix
# outstore$adjmatrix