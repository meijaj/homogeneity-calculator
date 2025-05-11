# Homogeneity calculations for CRMs 
# Author: Juris Meija, NRC Canada
# Date: 2019-2024, version Aug 2024

# Brief description
# This calculator fits the Bayesian normal-normal random effects model
# to the replicate measurement results of various Certified Reference Material (CRM) units

require(shiny)
require(shinyjs) # Javascript for shiny
require(shinyFeedback)
require(shinycssloaders)
require(rhandsontable) # Hot tables in shiny
require(rjags)
require(readxl) # read xlsx files
require(fitdistrplus) # fit gamma distribution
require(logspline) # smoothed kernel density

## TEST DATA
init.df = data.frame(unit=c('U1', 'U1', 'U1', 'U7', 'U7', 'U7', 'U23', 'U23'),
                     result=c(98, 99, 97, 95, 96, 94, 103, 102))

fit_gamma = function(m, k){
  require(rriskDistributions)
  get.gamma.par(p=c(0.10,0.50,0.90), q = m * c(1/k, 1, k), show.output = FALSE, plot = FALSE)
  
}

### server functions ####
server <- function(input, output, session) {
  
  v <- reactiveValues(# general parameters
                      dig = 5, failedconvergence=FALSE, init = TRUE,
                      # MCMC draws from the posterior distributions of u_meas and u_hom
                      MCMC_mu = 0, MCMC_meas=0, MCMC_hom=0, summary=NULL,
                      # data
                      unit = NULL, result = NULL, id = NULL, N = NULL, L = NULL,
                      # prior for mu
                      mu_prior = NULL, mu_prior_sd = NULL,
                      # prior for u_meas
                      u_meas_prior_median = NULL, u_meas_prior_gamma_shape= NULL, u_meas_prior_gamma_rate = NULL,
                      # prior for u_hom
                      u_hom_prior_median = NULL, u_hom_prior_gamma_shape= NULL, u_hom_prior_gamma_rate = NULL
                      )
  
  output$exampleData <- downloadHandler(
    filename = function() 'iso35-2017-appendix-c1.xlsx',
    content = function(file) {file.copy('iso35-2017-appendix-c1.xlsx', file) }
  )
  
  output$hot <- renderRHandsontable({
    infile <- input$upload
    
    if(is.null(infile)) { if(is.null(input$hot)) { DF = init.df } else { DF = hot_to_r(input$hot) } } 
    else {
      DF = read_xlsx(infile$datapath, col_names = TRUE, col_types = c('text','numeric'))[,1:2]
    }

    hide('f')
    hide('tablecaption')
    hide('table')
    hide('modeltext')
    
    names(DF) <- c('unit','result')
    rhandsontable(DF, readOnly = FALSE, stretchH = "all", selectCallback = TRUE) %>%
      hot_context_menu(allowColEdit = FALSE ) %>%
      hot_validate_numeric(cols=2) %>%
      hot_col(1:2, halign = 'htCenter')
  })
  
  observe({
    if (is.null(input$hot)) { z = init.df } else { z = hot_to_r(input$hot) }
    v$N_MCMC = input$N_MCMC
    dig = input$N_digits
    
    v$u_meas_pdf = 'half-Cauchy'
    v$u_hom_pdf = 'half-Cauchy'
    
    gp_meas = gp_hom = c(1, 1)
    
    v$unit = z[,1]
    v$id = as.factor(z[,1])
    v$result = as.double(z[,2])
    v$N = length(z[,1])
    v$L = length(unique(z[,1]))
    
    a = list()
    a$mu_prior = signif(median(sapply(split(z[,2], z[,1]), median, na.rm=TRUE), na.rm=TRUE), digits = dig)
    updateNumericInput(session, "mu_prior",    value = a$mu_prior )
    
    a$mu_prior_sd = signif(sd(z[,2], na.rm=TRUE), digits = dig)
    updateNumericInput(session, "mu_prior_sd", value = a$mu_prior_sd )
    
    a$u_meas_prior_median = signif(median(sapply(split(z[,2],  z[,1]), sd, na.rm=TRUE), na.rm=TRUE), digits =dig)
    updateNumericInput(session, "u_meas_prior_median", value = a$u_meas_prior_median )
    
    v$u_meas_prior_gamma_shape = gp_meas[1]
    v$u_meas_prior_gamma_rate  = gp_meas[2]
    
    a$u_hom_prior_median = signif(sd(sapply(split(z[,2],  z[,1]), median, na.rm=TRUE), na.rm=TRUE), digits =dig)
    updateNumericInput(session, "u_hom_prior_median", value = a$u_hom_prior_median )
    
    v$u_hom_prior_gamma_shape = gp_hom[1]
    v$u_hom_prior_gamma_rate  = gp_hom[2]
  })
  
  observeEvent(input$u_meas_pdf, {
    # estimate prior (gamma) distribution parameters
    gp_meas = c(1,1)
    v$u_meas_pdf = 'half-Cauchy'
    if(input$u_meas_pdf == 'Gamma'){
      v$u_meas_pdf = 'Gamma'
      gp_meas = fit_gamma(input$u_meas_prior_median, input$u_meas_prior_gamma_factor)
    }
    v$u_meas_prior_gamma_shape = gp_meas[1]
    v$u_meas_prior_gamma_rate  = gp_meas[2] 
  })
  
  observeEvent(input$u_hom_pdf, {
    # estimate prior (gamma) distribution parameters
    gp_hom = c(1,1)
    v$u_hom_pdf = 'half-Cauchy'
    if(input$u_hom_pdf == 'Gamma'){
      v$u_hom_pdf = 'Gamma'
      gp_hom = fit_gamma(input$u_hom_prior_median, input$u_hom_prior_gamma_factor)
    }
    v$u_hom_prior_gamma_shape = gp_hom[1]
    v$u_hom_prior_gamma_rate  = gp_hom[2]
  })
  
  v$init = FALSE
  
  observeEvent(input$N_digits, { if( input$N_digits %in% c(0:6) ) { hideFeedback("N_digits") } else { showFeedbackDanger("N_digits", "Invalid number of digits") }  })
  observeEvent(input$N_MCMC,   { if( input$N_MCMC >= 10000 ) { hideFeedback("N_MCMC") } else { showFeedbackDanger("N_MCMC", "Invalid number of MCMC draws") }  })
  observeEvent(input$mu_prior,  { 
    if( is.numeric(input$mu_prior) ) { 
      hideFeedback("mu_prior")
      v$mu_prior = input$mu_prior
      } else { showFeedbackDanger("mu_prior", "Invalid value") }
    })
  observeEvent(input$mu_prior_sd,     { 
    if( is.numeric(input$mu_prior_sd) & input$mu_prior_sd > 0 ) { 
      hideFeedback("mu_prior_sd") 
      v$mu_prior_sd = input$mu_prior_sd
      } else { showFeedbackDanger("mu_prior_sd", "Invalid value") }  
    })
  
  observeEvent(input$u_meas_prior_median, { 
    if( is.numeric(input$u_meas_prior_median) & as.double(input$u_meas_prior_median) > 0 ) { 
      hideFeedback("u_meas_prior_median")
      v$u_meas_prior_median = input$u_meas_prior_median
      
      gp_meas = c(1,1)
      v$u_meas_pdf = 'half-Cauchy'
      if(input$u_meas_pdf == 'Gamma'){
        v$u_meas_pdf = 'Gamma'
        gp_meas = fit_gamma(input$u_meas_prior_median, input$u_meas_prior_gamma_factor)
      }
      v$u_meas_prior_gamma_shape = gp_meas[1]
      v$u_meas_prior_gamma_rate  = gp_meas[2]
      
      } else { showFeedbackDanger("u_meas_prior_median", "Invalid value") }  
    })
  observeEvent(input$u_hom_prior_median,  { 
    if( is.numeric(input$u_hom_prior_median) & input$u_hom_prior_median > 0  )  { 
      hideFeedback("u_hom_prior_median") 
      v$u_hom_prior_median = input$u_hom_prior_median
      
      gp_hom = c(1,1)
      v$u_hom_pdf = 'half-Cauchy'
      if(input$u_hom_pdf == 'Gamma'){
        v$u_hom_pdf = 'Gamma'
        gp_hom = fit_gamma(input$u_hom_prior_median, input$u_hom_prior_gamma_factor)
      }
      v$u_hom_prior_gamma_shape = gp_hom[1]
      v$u_hom_prior_gamma_rate  = gp_hom[2]
      
      } else { showFeedbackDanger("u_hom_prior_median", "Invalid value") }  
    })
  
  observeEvent(input$u_meas_prior_gamma_factor,  { 
    if( is.numeric(input$u_meas_prior_gamma_factor) & input$u_meas_prior_gamma_factor >= 1.05  )  { 
      hideFeedback("u_meas_prior_gamma_factor") 
      v$u_meas_prior_gamma_factor = input$u_meas_prior_gamma_factor
      
      gp_meas = c(1,1)
      v$u_meas_pdf = 'half-Cauchy'
      if(input$u_meas_pdf == 'Gamma'){
        v$u_meas_pdf = 'Gamma'
        gp_meas = fit_gamma(input$u_meas_prior_median, input$u_meas_prior_gamma_factor)
      }
      v$u_meas_prior_gamma_shape = gp_meas[1]
      v$u_meas_prior_gamma_rate  = gp_meas[2]
      
    } else { showFeedbackDanger("u_meas_prior_gamma_factor", "Invalid value") }  
  })
  observeEvent(input$u_hom_prior_gamma_factor,  { 
    if( is.numeric(input$u_hom_prior_gamma_factor) & input$u_hom_prior_gamma_factor >= 1.05  )  { 
      hideFeedback("u_hom_prior_gamma_factor") 
      v$u_hom_prior_gamma_factor = input$u_hom_prior_gamma_factor
      
      gp_hom = c(1,1)
      v$u_hom_pdf = 'half-Cauchy'
      if(input$u_hom_pdf == 'Gamma'){
        v$u_hom_pdf = 'Gamma'
        gp_hom = fit_gamma(input$u_hom_prior_median, input$u_hom_prior_gamma_factor)
      }
      v$u_hom_prior_gamma_shape = gp_hom[1]
      v$u_hom_prior_gamma_rate  = gp_hom[2]
      
      
    } else { showFeedbackDanger("u_hom_prior_gamma_factor", "Invalid value") }  
  })
  
  observeEvent(input$button, {
    
    validate(
      need(input$N_digits %in% 0:6, 'DIGITS'),
      need(input$N_MCMC > 10000, 'N_MCMC'),
      need(is.numeric(input$mu_prior), 'MU_PRIOR'),
      need(is.numeric(input$mu_prior_sd) & input$mu_prior_sd > 0, 'U_MU_PRIOR'),
      need(is.numeric(input$u_meas_prior_median) & input$u_meas_prior_median > 0, 'U_MEAS_PRIOR_MEDIAN'),
      need(is.numeric(input$u_hom_prior_median) & input$u_hom_prior_median > 0, 'U_HOM_PRIOR_MEDIAN'),
      need(is.numeric(input$u_meas_prior_gamma_factor) & input$u_meas_prior_gamma_factor >= 1.05, 'U_MEAS_PRIOR_GAMMA_FACTOR'),
      need(is.numeric(input$u_hom_prior_gamma_factor) & input$u_hom_prior_gamma_factor >= 1.05, 'U_HOM_PRIOR_GAMMA_FACTOR'),
      need(length(unlist(do.call(rbind, input$hot$params$data)[,1])) == length(unlist(do.call(rbind, input$hot$params$data)[,2])), 'EMPTY DATA FIELDS')
    )
    
    umeas_model = c(paste0("meas ~ dt(0, 1.0/(", v$u_meas_prior_median, ")^2, 1)T(0,) \n"),
                    paste0("meas ~ dgamma(", round(v$u_meas_prior_gamma_shape,4),", ", round(v$u_meas_prior_gamma_rate,4),") \n")
                    )[v$u_meas_pdf == c('half-Cauchy','Gamma')]
    uhom_model = c(paste0("hom ~ dt(0, 1.0/(", v$u_hom_prior_median, ")^2, 1)T(0,) \n"),
                   paste0("hom ~ dgamma(", round(v$u_hom_prior_gamma_shape, 4),", ", round(v$u_hom_prior_gamma_rate,4),") \n")
                   )[v$u_hom_pdf == c('half-Cauchy','Gamma')]
    
    modelstring=paste("model {
      for(i in 1:N) { X[i] ~ dnorm(mu + theta[id[i]], 1.0/meas^2) }
      for(i in 1:L) { theta[i] ~ dnorm(0, 1.0/hom^2) } \n",
      paste0("mu ~ dt(", v$mu_prior, ",1.0/(", v$mu_prior_sd/sqrt(3), ")^2, 3) \n"),
      umeas_model, uhom_model, "} \n", collapse="\n")
    
    toggle(id = 'progress_text', condition = FALSE)
    output$progress_text = renderText('')
    show('progress')
    hide('f')
    hide('tablecaption')
    hide('table')
    hide('modeltext')
    
    model=jags.model(textConnection(modelstring), data=list(X=v$result, N=v$N, L = v$L, id=v$id), quiet=TRUE)
    update(model, n.iter=v$N_MCMC)
    output = coda.samples(model=model, variable.names=c("mu","meas", "hom"), n.iter = v$N_MCMC * 10, thin = 10)
    
    v$dig = input$N_digits
    
    # MCMC draws from the posterior pdfs
    v$MCMC_mu <- c(as.mcmc(output)[,'mu'])
    v$MCMC_meas <- c(as.mcmc(output)[,'meas'])
    v$MCMC_hom <- c(as.mcmc(output)[,'hom'])
    
    # summary output
    v$summary <- summary(output)
    
    # estimate gamma distribution that approximates the posterior distribution of u_hom
    v$hom_gamma <- fitdist(c(output[[1]][,'hom']), "gamma", method = "mle")$estimate
    
    # MCMC convergence test (Gewke)
    v$failedconvergence <- min(p.adjust(2*pnorm(-abs(geweke.diag(as.mcmc(output))[[1]])), method="BH")) < 0.05
    
    toggle(id = 'progress_text', condition = TRUE)
    hide('progress')
    show('f')
    show('tablecaption')
    show('table')
    show('modeltext')
  })
  
  output$f <- renderPlot({
    
    if(length(v$MCMC_meas)>1 & length(v$MCMC_hom)>1){  
    par(mfrow=c(1,2))
      qq <- quantile(v$MCMC_meas, c(0.01, 0.05, 0.95, 0.99, 0.50))
      dens <- density(v$MCMC_meas)
      
      # smoothed kernel density using splines to approximate the log-density
      m <- logspline(v$MCMC_meas, maxknots = 9)
      x0.logspline = seq(qq[1], qq[4], length.out = 300)
      x1.logspline = seq(qq[2], qq[3], length.out = 300)
      y0.logspline = sapply(x0.logspline, function(x) dlogspline(x, m))
      y1.logspline = sapply(x1.logspline, function(x) dlogspline(x, m))
      
      plot(m, main='Uncertainty due to measurement', pch='', col='steelblue3', 
           yaxt='n', ylab='', xlab='', xlim=c(qq[1], qq[4]))
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey98")
      # 5-95 % credible interval for the posterior
      polygon(c(x1.logspline, rev(x1.logspline)), c(y1.logspline, rep(0,length(x1.logspline))), col='#4F94CD7F', border='gray98')
      # prior
      xx <- seq(0, qq[4], length.out = 300)
      if(v$u_meas_pdf == 'half-Cauchy') lines(x=xx, y=dcauchy(xx, location=0, scale=v$u_meas_prior_median), col='tomato', lwd=4)
      if(v$u_meas_pdf == 'Gamma') lines(x=xx, y=dgamma(xx, shape=v$u_meas_prior_gamma_shape, rate=v$u_meas_prior_gamma_rate), col='tomato', lwd=4)
      # posterior
      lines(x=x0.logspline,y=y0.logspline, lwd=4, col='steelblue3')    
      
      segments(qq[2], 0, qq[3], 0)
      points(x=qq[5], y=0, pch=19, cex=1.3)
      legend("topright", legend = c('prior','posterior'), col = c('tomato','steelblue3'), lwd = c(4, 4), text.col = c('tomato','steelblue3'), text.font = c(2, 2), border = 'gray40')
      graphics::box()
      
      qq <- quantile(v$MCMC_hom, c(0.01, 0.05, 0.95, 0.99, 0.50))
      dens <- density(v$MCMC_hom)
      
      # smoothed kernel density using splines to approximate the log-density
      m <- logspline(v$MCMC_hom, maxknots = 9)
      x0.logspline = seq(qq[1],qq[4],length.out = 300)
      x1.logspline = seq(qq[2],qq[3],length.out = 300)
      y0.logspline = sapply(x0.logspline, function(x) dlogspline(x, m))
      y1.logspline = sapply(x1.logspline, function(x) dlogspline(x, m))
      
      plot(m, main='Uncertainty due to homogeneity', pch='', col='steelblue3', yaxt='n', ylab='', xlab='', xlim=c(qq[1], qq[4]))
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey98")
      # 5-95 % credible interval for the posterior
      polygon(c(x1.logspline, rev(x1.logspline)), c(y1.logspline, rep(0,length(x1.logspline))), col='#4F94CD7F', border='gray98')
      # prior
      xx <- seq(0, qq[4], length.out = 300)
      if(v$u_hom_pdf == 'half-Cauchy') lines(x=xx, y=dcauchy(xx, location=0, scale=v$u_hom_prior_median), col='tomato', lwd=4)
      if(v$u_hom_pdf == 'Gamma') lines(x=xx, y=dgamma(xx, shape=v$u_hom_prior_gamma_shape, rate=v$u_hom_prior_gamma_rate), col='tomato', lwd=4)
      
      # posterior
      lines(x=x0.logspline, y=y0.logspline, lwd=4, col='steelblue3')    
      
      segments(qq[2], 0, qq[3], 0)
      points(x=qq[5], y=0, pch=19, cex=1.3)
      legend("topright", legend = c('prior','posterior'), col = c('tomato','steelblue3'), lwd = c(4, 4), text.col = c('tomato','steelblue3'), text.font = c(2, 2), border = 'gray40')
      graphics::box()
      
    }
    })
 
  modeltext <- eventReactive(input$button, { 
    HTML(paste('<b>Statistical model</b>',
               'This calculator performs the Bayesian analysis of variance (BANOVA) using the following statistical measurement model, written in the BUGS programming language:',
               paste0('&nbsp;&nbsp;&nbsp; result[<i>i</i>] ~ normal(mean = <i>m</i> + <i>h</i>[unit[i]], sd = <i>u</i><sub>meas</sub>) for each observation <i>i</i> (<i>i</i> = 1...', v$N,")"),
               paste0('&nbsp;&nbsp;&nbsp; <i>h</i>[<i>u</i>] ~ normal(mean = 0, sd = <i>u</i><sub>hom</sub>) for each CRM unit <i>u</i> (<i>u</i> = 1...', v$L,")"),
               'The following priors were used for the model parameters:',
               paste0("&nbsp;&nbsp;&nbsp; <i>m</i> ~ t(mean = ", v$mu_prior, ", sd = ", v$mu_prior_sd,", df = 3)"),
               paste(ifelse(v$u_meas_pdf=="half-Cauchy",
                      paste0("&nbsp;&nbsp;&nbsp; <i>u</i><sub>meas</sub> ~ halfCauchy(median = ", v$u_meas_prior_median,")"),
                      paste0("&nbsp;&nbsp;&nbsp; <i>u</i><sub>meas</sub> ~ Gamma(shape = ", signif(v$u_meas_prior_gamma_shape, v$dig),", rate = ", signif(v$u_meas_prior_gamma_rate, v$dig), ")")
                      ),""),
               paste0(ifelse(v$u_hom_pdf=="half-Cauchy",
                      paste0("&nbsp;&nbsp;&nbsp; <i>u</i><sub>hom</sub> ~ halfCauchy(median = ", v$u_hom_prior_median,")"),
                      paste0("&nbsp;&nbsp;&nbsp; <i>u</i><sub>hom</sub> ~ Gamma(shape = ", signif(v$u_hom_prior_gamma_shape, v$dig),", rate = ", signif(v$u_hom_prior_gamma_rate, v$dig), ")")
               ),""),
               paste0("The model was fit to data using Markov-chain Monte Carlo method in R using rjags"),
               paste0("The (posterior) probability distribution for <i>u</i><sub>hom</sub> can be summarized using gamma distribution with shape = ",formatC(v$hom_gamma[1],digits=input$N_digits, format='f')," and rate = ", formatC(v$hom_gamma[2],digits=input$N_digits, format='f')),
               paste0(ifelse(v$failedconvergence,"<p style='color:red'>WARNING: MCMC may not have reached equilibrium. Results are not reliable. Try increasing the number of iterations.</p>","")),
               sep="<br/>"))
    })
  
  tablecaption <- eventReactive(input$button, {
    HTML(paste('<b>Summary of the results</b>'))
  })
  
  output$modeltext <- renderUI({  modeltext() })
  output$tablecaption <- renderUI({  tablecaption() })
  output$table <- renderTable( if(input$button & !is.null(v$MCMC_meas) & !is.null(v$MCMC_hom)){  
    
    tt = rbind(  quantile(v$MCMC_mu, c(0.05, 0.50, 0.95) ),
                 quantile(v$MCMC_meas, c(0.05, 0.50, 0.95) ),
                 quantile(v$MCMC_hom, c(0.05, 0.50, 0.95) ) )
    rownames(tt) = c('mu', 'u_meas', 'u_hom')
    colnames(tt) = c('5% quantile', 'median', '95% quantile')
    tt
  }, rownames = TRUE, digits = isolate(v$dig-1) ) 
  shinyjs::onclick("toggleextra", shinyjs::toggle(id = "filterextra", anim = TRUE))
  
  }

### user interface ####
ui <- fluidPage(
  useShinyFeedback(),
  titlePanel( title="NRC CRM Homogeneity uncertainty calculator" ),
  shinyjs::useShinyjs(),
  sidebarLayout(
    
    sidebarPanel(
      h5(tags$b("Enter (paste) the observed results")),
      helpText("right-click to add or delete rows"),textOutput("error"),
      rHandsontableOutput("hot", height = 300),
      h5(tags$b("Upload dataset")),
      helpText("The first two named columns of the data file will be taken for 'unit' and 'result', respectively."),
      fileInput("upload", label=NULL, placeholder="Upload data file (xlsx)", multiple = FALSE, accept = ".xlsx"),
      downloadLink("exampleData", label = "Download ISO 35:2017 example dataset (Appendix C.1)"),
      br(),
      h5(tags$b("Additional model settings "), a(id = "toggleextra", "show/hide")),
      shinyjs::hidden(div(id = "filterextra",
                          fluidRow(
                            column(6, numericInput("N_MCMC",  label = "Number of MCMC draws", value = 100000, min = 10000, step = 10000, width =  '85%')),
                            column(6, numericInput("N_digits", label = "Decimal digits for display", value = 3, max = 6, min = 0, step = 1, width='85%')),
                            column(12, h5('Information about the prior distributions of model parameters')),
                            column(6, numericInput("mu_prior",    label = "Student-t3 mean for m", value = 1, width = '85%')),
                            column(6, numericInput("mu_prior_sd", label = "Student-t3 sd for m", value = 100, min=0, width = '85%')),
                            column(6, selectInput("u_meas_pdf", label = "Select PDF for u_meas", choices = c("half-Cauchy", "Gamma"), width= '85%')),
                            column(6, selectInput("u_hom_pdf",  label = "Select PDF for u_hom", choices = c("half-Cauchy", "Gamma"), width= '85%')),
                            column(6, numericInput("u_meas_prior_median", label = "median for u_meas", value = 1, min=0, width = '85%')),
                            column(6, numericInput("u_hom_prior_median",  label = "median for u_hom", value = 1, min=0, width = '85%')),
                            column(12, h6('Type a sufficiently large value if u_hom or u_meas are trully unknown.')),
                            column(6, numericInput("u_meas_prior_gamma_factor", label = "reliability factor for u_meas", value = 9,  min = 1.05, max=100, width = '85%')),
                            column(6, numericInput("u_hom_prior_gamma_factor",  label = "reliability factor for u_hom", value = 100, min = 1.05, max=100, width = '85%')),
                            column(12, h6('When gamma distribution is chosen, approx. 80% of values for either u_hom or u_meas will be in the interval from median * 1/k to median * k. The reliability factors are ignord when half-Cauchy distribution is chosen.'))
                          ))),
      br(),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                       actionButton("button", label = "Fit model", icon = icon('bar-chart-o'))),
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",   
                       actionButton("button", label = "busy...", icon = icon('hourglass-half'))),
      br(),
      p("Juris Meija (2024) NRC CRM homogeneity calculator, v3")
    ),
    
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      br(),
      hidden(div(id = 'progress', withSpinner(textOutput(outputId = "progress_text")))),
      plotOutput("f"),
      htmlOutput("tablecaption"),
      tableOutput("table"),
      htmlOutput("modeltext"),br()
    )
  )
)

shinyApp(ui = ui, server = server)