############
## new class extending the ExpressionSet from Bioconductor!
## just containing one additional matrix with aggregated residual information per probe set.
setClass("ResExpressionSet",contains="ExpressionSet",representation( residuals="matrix" ),      prototype=list(residuals=matrix())
)

###############################
# new
if( !isGeneric("residuals") )
	setGeneric("residuals", function( object,...)
	standardGeneric("residuals"))


##################################################################################################
setMethod("residuals","ResExpressionSet",
          function( object, ... ){
              return( object@residuals )
          }
          )

newResExpressionSet <- function( x ){
    if( class( x )=="ExpressionSet" ){
        x <- as( x, "ResExpressionSet" )
        Res <- matrix( ncol=ncol( exprs( x ) ), nrow=nrow( exprs( x ) ) )
        colnames( Res ) <- colnames( exprs( x ) )
        rownames( Res ) <- rownames( exprs( x ) )
        x@residuals <- Res
        return( x )
    }else{
        stop( "Not yet implemented" )
    }
}
