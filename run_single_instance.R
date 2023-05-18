run_single_instance = function(file.date,i_day,poll)
{

##################################
## Using the SPDE approach and the INLA algorithm for spatio-temporal
## modelling and mapping
##
## adapted from M.Cameletti, F.Lindgren, D.Simpson, H.Rue
##################################

## User-overrideable configuration:
## Set any of these variables to alternate values before sourcing this file.

## To force a rerun of inla, set force.rerun=TRUE
## If "result" does not exist, inla will be rerun regardless of this setting.
if (!exists("force.rerun")) { force.rerun=TRUE }
## Do you want to run inla() on a remote server?
## Requires having run  inla.remote()  previously.
if (!exists("do.remote")) { do.remote = FALSE }
## Do you want to print (to png/pdf files)?
#if (!exists("do.print")) { do.print = FALSE }
do.print = FALSE
## Path for the data directory
if (!exists("data.path")) { data.path = "INLA_Data/" }

options(error=recover)

graphics.off()
if (force.rerun || !exists("result"))
{

    ## ################################
    ## Load the data
    ## ################################
    ##--- for the Nest stations and P days
    filename <- paste(data.path,"Italy_", poll, "data_", file.date, "_est.csv",sep="")
    OK <- file.exists(filename)
    if ( !(OK) ) {
       return(paste("File", filename, "not found",sep=" "))
    }
    print(paste("Using data from ", file.date), sep="")
    Italy_data_est <-read.table(paste(data.path,"Italy_", poll, "data_", file.date, "_est.csv",sep=""),header=TRUE,sep=",")
    Italy_data_val <-read.table(paste(data.path,"Italy_", poll, "data_", file.date, "_val.csv",sep=""),header=TRUE,sep=",")
    coords_est <-read.table(paste(data.path,"coords_est.csv",sep=""),header=TRUE,sep=",")
    coords_val <-read.table(paste(data.path,"coords_val.csv",sep=""),header=TRUE,sep=",")

    rownames(coords_est) = coords_est[,"Station.ID"]
    rownames(coords_val) = coords_val[,"Station.ID"]

    which_date = unique(Italy_data_est$Date)[i_day]
    print(paste("**---- You will get a prediction for ", which_date, "---**"))

    ## ################################
    ## Work out how many stations and days there are
    ## ################################
    n_stations_est <- length(coords_est$Station.ID)
    n_stations_val <- length(coords_val$Station.ID)
    n_data <- length(Italy_data_est$Station.ID)
    n_days <- as.integer(n_data/n_stations_est)
    n_cols <- ncol(Italy_data_est)

    ##--- Check that the data is OK
    if (n_data %% n_stations_est != 0) {
print(n_data)
print(n_stations_est)
        print("The number of data points needs to be an integer multiple of the number of stations!")
        return
    }

    ##--- Check validation data is OK
    if (length(Italy_data_val$Station.ID) != n_stations_val*n_days) {
        print("Something is wrong with the dimension of the validation data!")
        return
    }

    cat("# of stations used for estimation:", n_stations_est, "\n")
    cat("# of stations used for validation:", n_stations_val, "\n")
    cat("# of days:", n_days, "\n")
    cat("# of covariates:", n_cols-3, "\n")
    cat("list of covariates: ")
    for (i in 4:n_cols) {
     cat(colnames(Italy_data_est[i])," ")
    }
    cat("\n")

    ## ################################
    ## standardize covariates and calculate log of data
    ## ################################
    Italy_data_est$ENS = log(Italy_data_est$ENS)
    Italy_data_val$ENS = log(Italy_data_val$ENS)

    mean_covars = apply(Italy_data_est[,4:n_cols],2,mean)
    sd_covars = apply(Italy_data_est[,4:n_cols],2,sd)

    Italy_data_est[,4:n_cols] = scale(Italy_data_est[,4:n_cols], mean_covars, sd_covars)
    Italy_data_val[,4:n_cols] = scale(Italy_data_val[,4:n_cols], mean_covars, sd_covars)

    Italy_data_est$logData = log(Italy_data_est$Data)
    Italy_data_val$logData = log(Italy_data_val$Data)

    Italy_data_est$time = rep(1:n_days,each = n_stations_est)
    Italy_data_val$time = rep(1:n_days,each = n_stations_val)

    ## ################################
    ## Estimation
    ## ################################


    filename <- paste(data.path,"Bounds.Rdata",sep="")
    OK <- file.exists(filename)
    if ( OK ) {
       cat("Load previously saved triangulation ...\n")
       load(filename)
    }
    else {
       bounds.exist <- exists("bounds")
       if ( ! bounds.exist )
       {
         cat("Estimate triangulation ...\n")
         ## ################################
         ## Triangulation using coords_est
         ## ################################
         bounds <- inla.nonconvex.hull(cbind(coords_est$UTMX,coords_est$UTMY), convex = 55.0)
         mesh = inla.mesh.2d(boundary = bounds, cutoff = 10.0, max.edge = c(40, 120), offset=cbind(300,500))
       }
       save(bounds, mesh, file=filename)
    }

    ##--- Plot the triangulation
    if ( do.print )
    {
        ##--- Borders of Italy (in km)
        # Sicily
        borders1 <-read.table(paste(data.path,"Italian_boundaries1.csv",sep=""),header=TRUE,sep=",")
        # Italian peninsula
        borders2 <-read.table(paste(data.path,"Italian_boundaries2.csv",sep=""),header=TRUE,sep=",")
        # Sardinia
        borders3 <-read.table(paste(data.path,"Italian_boundaries3.csv",sep=""),header=TRUE,sep=",")

        inla.dev.new()
        pdf("triangulation_italy.pdf")
        cat("plot triangulation ...\n")
        plot(mesh)
        lines(borders1, lwd=1.5, col="black")
        lines(borders2, lwd=1.5, col="black")
        lines(borders3, lwd=1.5, col="black")
        points(coords_est$UTMX, coords_est$UTMY, pch=20, cex=0.5, col="blue")
        points(coords_val$UTMX, coords_val$UTMY, pch=20, cex=0.5, col="red")
    }
    #if (do.print) { dev.off() } else { inla.dev.new() }

    ## ################################
    ## Make the SPDE object and the formula
    ## ################################

    ##--- Construct the SPDE object
    spde = inla.spde2.matern(mesh=mesh, alpha=2)

    ##--- Observation structure for estimation data
    A.est =
        inla.spde.make.A(mesh,
                         loc=
                         as.matrix(coords_est[Italy_data_est$Station.ID, c("UTMX","UTMY")]),
                         group=Italy_data_est$time,
                         n.group=n_days
                         )
    ##--- Observation structure for validation data
    A.val =
        inla.spde.make.A(mesh,
                         loc=
                         as.matrix(coords_val[Italy_data_val$Station.ID, c("UTMX","UTMY")]),
                         group=Italy_data_val$time,
                         n.group=n_days
                         )

    ##--- Observation structure for field prediction
    A.pred =
        inla.spde.make.A(mesh, group=i_day, n.group=n_days)

    field.indices =
        inla.spde.make.index("field",
                             n.spde=spde$n.spde,
                             n.group=n_days)

    stack.est =
        inla.stack(data=list(logData=Italy_data_est$logData),
                   A=list(A.est, 1),
                   effects=
                   list(c(field.indices, list(Intercept=1)), list(Italy_data_est[,4:n_cols])), tag="est")

    stack.val =
        inla.stack(data=list(logData=NA),
                   A=list(A.val, 1),
                   effects=
                   list(c(field.indices, list(Intercept=1)), list(Italy_data_val[,4:n_cols])), tag="val")

    scaled.mesh.loc =
        list(UTMX=(rep(scale(mesh$loc[,1], mean_covars["UTMX"], sd_covars["UTMX"]), n_days)),
             UTMY=(rep(scale(mesh$loc[,2], mean_covars["UTMY"], sd_covars["UTMY"]), n_days)))

    stack.pred =
        inla.stack(data=list(logData=NA),
                   A=list(A.pred),
                   effects=
                   list(c(field.indices, scaled.mesh.loc, list(Intercept=1))), tag="pred")

    stack = inla.stack(stack.est, stack.val, stack.pred)

    formula <- (logData ~ -1 + Intercept + ENS + WS + T2 + PBL00 + PBL12
                             + TP + POP + ALT + Imperviousness + Builtup + roads + highways
                             + HUD + LUD + Agriculture + Forests
              + f(field, model=spde, group=field.group, control.group=list(model="ar1"))
               )
    print(formula)

    ## ################################
    ## Call INLA and get results
    ## ################################
    ## Attention!!! Set do.remote=FALSE to disable inla.remote
    ## reordering="metis" is needed to prevent crashes on some systems.
    filename <- paste(data.path,"INLA_",poll,"_",file.date,".Rdata",sep="")
    OK <- file.exists(filename)
    if ( OK ) {
       cat("Load previously saved INLA data ...\n")
       load(filename)
    }
    else { 
      cat("call INLA ...\n")
      if (!exists("do.remote") || do.remote) {
          result =
              inla(formula,
                   data=inla.stack.data(stack, spde=spde),
                   family="gaussian",
                   control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
                   control.compute=list(cpo=FALSE),
                   control.inla = list(reordering = "metis"),
                   keep=FALSE, verbose=TRUE,
                   inla.call="remote", num.threads=12)
      } else {
          result =
              inla(formula,
                   data=inla.stack.data(stack, spde=spde),
                   family="gaussian",
                   control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
                   control.compute=list(cpo=FALSE),
                   control.inla = list(reordering = "metis"),
                   keep=FALSE, verbose=FALSE)
      }
      save(result, file=filename)      
    }
    print(summary(result))

} else {
    print.noquote("Note: Detected existing 'result' variable,")
    print.noquote("      and force.rerun=FALSE, so not rerunning inla.")
    print.noquote("      Set force.rerun=TRUE to recompute the results.")
}

## ############################
## Extract results
## ############################
cat("extract results ...\n")

##--- Results for fixed effects - covariate coeffs -
beta = result$summary.fixed[,"mean"]
beta_sd = result$summary.fixed[,"sd"]

## Extract data information from estimation data
estimation = list()
index = inla.stack.index(stack,"est")$data
tmp.mean <- result$summary.linear.predictor[index,"mean"]
tmp.sd   <- result$summary.linear.predictor[index,"sd"]
estimation$res = exp(Italy_data_est$logData) - exp(tmp.mean)
Q0 <- qlnorm(0.025, meanlog=tmp.mean, sdlog=tmp.sd)
Q1 <- qlnorm(0.250, meanlog=tmp.mean, sdlog=tmp.sd)
Q3 <- qlnorm(0.750, meanlog=tmp.mean, sdlog=tmp.sd)
Q4 <- qlnorm(0.975, meanlog=tmp.mean, sdlog=tmp.sd)

#RMSE
est_res <- sqrt(mean((Italy_data_est[["Data"]]-exp(Italy_data_est[["ENS"]]*sd_covars[1]+mean_covars[1]))^2,na.rm=T)) # raw ensemble data
print(paste("average residual for estimation data before optimization", est_res))
est_res <- sqrt(mean((estimation$res)^2,na.rm=T)) #spde
print(paste("average residual for estimation data after optimization ", est_res))
#CORR
corrcoef1 <- cor(Italy_data_est["Data"], exp(Italy_data_est["ENS"]*sd_covars[1]+mean_covars[1]), use="pairwise.complete.obs", method="pearson")
print(paste("corr coeff for estimation data before optimization", corrcoef1))
corrcoef2 <- cor(exp(Italy_data_est$logData), exp(tmp.mean), use="pairwise.complete.obs", method="pearson")
print(paste("corr coeff for estimation data after optimization ", corrcoef2))

# export table
ExportedResults <- matrix(nrow=length(Italy_data_est$logData),ncol=9)
ExportedResults[,1] <- Italy_data_est$Date
ExportedResults[,2] <- Italy_data_est$Station.ID
ExportedResults[,3] <- Italy_data_est$Data
ExportedResults[,4] <- exp(Italy_data_est[["ENS"]]*sd_covars[1] + mean_covars[1])
ExportedResults[,5] <- exp(tmp.mean)
ExportedResults[,6] <- Q0
ExportedResults[,7] <- Q1
ExportedResults[,8] <- Q3
ExportedResults[,9] <- Q4
colnames(ExportedResults) <- c("Date", "Station.ID", "Data", "ENS", "Estimated", "Q0", "Q1", "Q3", "Q4")
ExportedRes1 <- ExportedResults[ExportedResults[,"Date"]!=which_date]
ExportedRes2 <- ExportedResults[ExportedResults[,"Date"]==which_date]
ncol_ExpRes <- 9
nrow_ExpRes <- NROW(ExportedRes1)/ncol_ExpRes
dim(ExportedRes1) <- c(nrow_ExpRes, ncol_ExpRes)
colnames(ExportedRes1) <- c("Date", "Station.ID", "Data", "ENS", "Estimated", "Q0", "Q1", "Q3", "Q4")
fileout=paste(data.path,"Italy_", poll, "data_", file.date, "_EstimationResults.csv",sep="")
write.table(ExportedRes1, fileout, append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

nrow_ExpRes <- NROW(ExportedRes2)/ncol_ExpRes
dim(ExportedRes2) <- c(nrow_ExpRes, ncol_ExpRes)
colnames(ExportedRes2) <- c("Date", "Station.ID", "Data", "ENS", "Estimated", "Q0", "Q1", "Q3", "Q4")
fileout=paste(data.path,"Italy_", poll, "data_", file.date, "_PredictedResults.csv",sep="")
write.table(ExportedRes2, fileout, append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)


## Extract information validation data
validation = list()
index = inla.stack.index(stack,"val")$data
tmp.mean = result$summary.linear.predictor[index,"mean"]
tmp.sd = result$summary.linear.predictor[index,"sd"]
validation$res = exp(Italy_data_val$logData) - exp(tmp.mean)
validation$logres = Italy_data_val$logData - tmp.mean
validation$res.std = validation$logres /
    sqrt(tmp.sd^2 + 1/result$summary.hyperpar[1,"mean"])
validation$p = pnorm(validation$res.std)
validation$rmse = sqrt(mean(validation$logres^2, na.rm=TRUE))
validation$cor = cor(Italy_data_val$logData, tmp.mean,
                     use="pairwise.complete.obs", method="pearson")
validation$cover = mean((validation$p>0.025)&(validation$p<0.975), na.rm=TRUE)

Q0 <- qlnorm(0.025, meanlog=tmp.mean, sdlog=tmp.sd)
Q1 <- qlnorm(0.250, meanlog=tmp.mean, sdlog=tmp.sd)
Q3 <- qlnorm(0.750, meanlog=tmp.mean, sdlog=tmp.sd)
Q4 <- qlnorm(0.975, meanlog=tmp.mean, sdlog=tmp.sd)

## write validation results on output file
#file.out=paste(data.path, "validation_", file.date, ".csv",sep="")
#write.matrix(paste(exp(Italy_data_val$logData), exp(tmp.mean), sep=","),file=file.out)

#RMSE
val_res <- sqrt(mean((Italy_data_val[["Data"]]-exp(Italy_data_val[["ENS"]]*sd_covars[1]+mean_covars[1]))^2,na.rm=T)) # raw ensemble data
print(paste("average residual for validation data before optimization", val_res))
val_res <- sqrt(mean((validation$res)^2,na.rm=T)) #spde
print(paste("average residual for validation data after optimization ", val_res))
#CORR
corrcoef1 <- cor(Italy_data_val["Data"], exp(Italy_data_val["ENS"]*sd_covars[1]+mean_covars[1]), use="pairwise.complete.obs", method="pearson")
print(paste("corr coeff for validation data before optimization", corrcoef1))
corrcoef2 <- cor(exp(Italy_data_val$logData), exp(tmp.mean), use="pairwise.complete.obs", method="pearson")
print(paste("corr coeff for validation data after optimization ", corrcoef2))

# export table
ExportedResults <- matrix(nrow=length(Italy_data_val$logData),ncol=9)
ExportedResults[,1] <- Italy_data_val$Date
ExportedResults[,2] <- Italy_data_val$Station.ID
ExportedResults[,3] <- Italy_data_val$Data
ExportedResults[,4] <- exp(Italy_data_val[["ENS"]]*sd_covars[1] + mean_covars[1])
ExportedResults[,5] <- exp(tmp.mean)
ExportedResults[,6] <- Q0
ExportedResults[,7] <- Q1
ExportedResults[,8] <- Q3
ExportedResults[,9] <- Q4
colnames(ExportedResults) <- c("Date", "Station.ID", "Data", "ENS", "Estimated", "Q0", "Q1", "Q3", "Q4")
ExportedResults <- ExportedResults[ExportedResults[,"Date"]!=which_date]
ncol_ExpRes <- 9
nrow_ExpRes <- NROW(ExportedResults)/ncol_ExpRes
dim(ExportedResults) <- c(nrow_ExpRes, ncol_ExpRes)
colnames(ExportedResults) <- c("Date", "Station.ID", "Data", "ENS", "Estimated", "Q0", "Q1", "Q3", "Q4")
fileout=paste(data.path,"Italy_", poll, "data_", file.date, "_ValidationResults.csv",sep="")
write.table(ExportedResults, fileout, append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

#plots
if ( do.print )
{
  inla.dev.new()
  par(mfrow=c(3,1))
  plot(Italy_data_val$logData,t="l",ylim=c(1,5.5))
  plot(tmp.mean, col=2,t="l",ylim=c(1,5.5))
  inla.dev.new()
  par(mfrow=c(1,1))
  plot(Italy_data_val$logData,tmp.mean, t="p",xlim=c(1,5.5),ylim=c(1,5.5))
  inla.dev.new()
  par(mfrow=c(1,1))
  plot(exp(Italy_data_val$logData),exp(tmp.mean), t="p",xlim=c(0,100),ylim=c(0,100))
}

##--- Extract SPDE results
result.field = inla.spde.result(result, "field", spde, do.transform=TRUE)


field_mean = matrix(result.field$summary.values$mean, mesh$n, n_days)
field_sd   = matrix(result.field$summary.values$sd, mesh$n, n_days)
field_pred_mean =
    result$summary.linear.predictor[inla.stack.index(stack,"pred")$data, "mean"]
field_pred_sd =
    result$summary.linear.predictor[inla.stack.index(stack,"pred")$data, "sd"]

##--- ar1 parameter (if any)
result$summary.hyperpar["GroupRho for field",]

##--- sigma2eps (1/precision)
sigma2eps_marg =
    inla.tmarginal(function(x) 1/x,
                   result$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2eps_m1 = inla.emarginal(function(x) x, sigma2eps_marg)
sigma2eps_m2 = inla.emarginal(function(x) x^2, sigma2eps_marg)
sigma2eps_stdev = sqrt(sigma2eps_m2 - sigma2eps_m1^2)
sigma2eps_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2eps_marg)

var.nom.marg = result.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)

range.nom.marg = result.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)

approx.hyperpar =
    rbind(obs.var=c(sigma2eps_m1, sigma2eps_stdev, sigma2eps_quantiles),
          spde.var.nom=c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
          spde.range.nom=c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
          AR.rho=result$summary.hyperpar["GroupRho for field",1:5])
print(approx.hyperpar)


##--- Load the full grid and covars
Griglia <- read.table(paste(data.path,"Griglia.csv",sep=""),sep=",",header=TRUE)

filename <- paste(data.path,"Italy_", poll, "data_", file.date, "_grid.csv",sep="")
OK <- file.exists(filename)
if ( !(OK) ) {
  return(paste("File", filename, "not found",sep=" "))
}
covars <- read.table(filename,header=TRUE,sep=",")

n_cols_covars <- ncol(covars)
n_rows_covars <- nrow(covars)

# log and normalize covars
covars[,3] <- log(covars[,3])
covars[,3:n_cols_covars] = scale(covars[,3:n_cols_covars], mean_covars, sd_covars)

##--- Create mapping between mesh and grid
proj_grid =     
    inla.mesh.projector(mesh,
                        xlim=range(Griglia[,2]),
                        ylim=range(Griglia[,3]),
                        dims=c(384,339))

##--- Project the latent field + geographical covariates
grid_latent_mean = inla.mesh.project(proj_grid, field_pred_mean)
grid_latent_sd = inla.mesh.project(proj_grid, field_pred_sd)

grid_mean <- matrix(0.0, n_rows_covars, 1)
grid_var <- matrix(0.0, n_rows_covars, 1)
prediction_points <- Griglia[,1]

##--- Reconstruct the logData field (trend + latent field)
##--- Assuming independence between all non-geographical components
k <- 0
ExportedResults <- matrix(nrow=n_rows_covars,ncol=12)
for (n in c(1:length(prediction_points))) {
  if ( !(is.na(prediction_points[n])) ) {
   k <- k + 1
   grid_mean[k] <- grid_latent_mean[n]
   grid_var[k] <- grid_latent_sd[n]
   grid_mean[k] <- grid_mean[k] + covars[k,3]*beta[2]    # ENS
   grid_var[k] <- grid_var[k] + (covars[k,3]*beta[2])^2
   for (i in c(3:length(beta)) ) {
    grid_mean[k] <- grid_mean[k] + covars[k,i+3]*beta[i]
    grid_var[k] <- grid_var[k] + (covars[k,i+3]*beta[i])^2
   }
   ExportedResults[k,1] <- n-1
   ExportedResults[k,2] <- 1000.0*Griglia[n,2]
   ExportedResults[k,3] <- 1000.0*Griglia[n,3]
   ExportedResults[k,4] <- exp(grid_mean[k])
   ExportedResults[k,5] <- exp(covars[k,3]*sd_covars[1] + mean_covars[1])
  }
}
grid_sd <- sqrt(grid_var)

ExportedResults[,6]  <- qlnorm(0.025, meanlog=log(ExportedResults[,4]), sdlog=log(ExportedResults[,5]))
ExportedResults[,7]  <- qlnorm(0.250, meanlog=log(ExportedResults[,4]), sdlog=log(ExportedResults[,5]))
ExportedResults[,8]  <- qlnorm(0.500, meanlog=log(ExportedResults[,4]), sdlog=log(ExportedResults[,5]))
ExportedResults[,9]  <- qlnorm(0.750, meanlog=log(ExportedResults[,4]), sdlog=log(ExportedResults[,5]))
ExportedResults[,10] <- qlnorm(0.975, meanlog=log(ExportedResults[,4]), sdlog=log(ExportedResults[,5]))

u_level <- log(40.0)
ExportedResults[,11] <- pnorm((grid_mean-u_level)/grid_sd)
u_level <- log(50.0)
ExportedResults[,12] <- pnorm((grid_mean-u_level)/grid_sd)


# export table
colnames(ExportedResults) <- c("id", "UTMX", "UTMY", poll, "ENS", "Q0", "Q1", "QM", "Q3", "Q4", "P40", "P50")
fileout=paste(data.path,"Italy_", poll, "data_", file.date, "_PredictionGrid.csv",sep="")
write.table(ExportedResults, fileout, append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

}
