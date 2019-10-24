# Ejemplo de identificación, estimación,  diagnósticos y pronósticos en un modelo ARIMA
#  
# ----------------------------------------------------------------------------
rm(list=ls(all=TRUE))  # elimina todos los objetos creados que todavía estén en R

# Cargar librerías
library(fpp2)
library(TSA)
library(car)
library(FitAR)
library(urca)            
library(fBasics)

# Entrada de datos desde un archivo de texto
(z=ts(scan("D:/Curso_Econometria_Avanz_EAFIT/datos/PIB_real_2000_1_2015_4.txt")))

# ETAPA DE IDENTIFICACIÓN
# gráfica de la serie
plot.ts(z, type="o", cex=0.6)
#
# Análisis de estabilidad de la varianza: Es necesario usar la
# transformación Box-Cox?

# Análisis gráfico
par(mfrow=c(2,1))
plot.ts(z, type="o", cex=0.6)
plot.ts(diff(z), type="o", cex=0.6)

# usando la librería forecast: supone que los datos son dependientes
BoxCox.lambda(z, method = c("loglik"), lower = -2, upper = 2)  # maximiza log likelihood
# sugiere una transformación logaritmica

par(mfrow=c(2,1))
plot.ts(log(z), type="o", cex=0.6)
plot.ts(diff(log(z)), type="o", cex=0.6)

# otra alernativa usando la librería forecast
BoxCox.lambda(z, method = c("guerrero"), lower = -2, upper = 2)  #minimiza coef. de variación

# usando la librería FitAR: supone que los datos son dependientes
(p=round(length(z)^(1/3)))
mod=arima(z, c(p, 1, 0))
tsdiag(mod)                 # algunos diagnósticos 
BoxCox(mod)

# transformación: lambda=1, es decir ninguna  
# se sigue trabajando con la serie z sin transformación Box-Cox

# determinación del valor de d
# correlogramas muestrales de la serie original
(max_rezag=round(length(z)/4))
par(mfrow=c(2,1))
Acf(z, lag.max=max_rezag, ci=0,ylim=c(-1,1))    # ACF decrece lentamente
Pacf(z, lag.max=max_rezag, ci=0, ylim=c(-1,1))  # PACF concentrada en el primer rezago

# La serie es no estacionaria: Qué tipo de serie es: TS o DS?
# prueba de raíces unitarias
# usando la librería urca
(maxlag=floor(12*(length(z)/100)^(1/4)))
ru_z=ur.df(z, type = c("trend"), lags=maxlag, selectlags = c("AIC")) # Usando AIC
summary(ru_z)
# valide el modelo

ru_z=ur.df(z, type = c("trend"), lags=maxlag, selectlags = c("BIC")) # Usando BIC
summary(ru_z)
# valide el modelo

# redefinición el modelo basada en el criterio BIC
ru_z=ur.df(z, type = c("trend"), lags=0)   # redefinición del modelo
summary(ru_z)
# valide el modelo

# se chequea cada modelo y se usa la prueba usando el modelo sugerido por el criterio AIC
  # validación de la ecuación ADF
    resid=ru_z@testreg$residuals             # residuales del modelo ADF
    plot(ru_z)                               # diagnósticos 
    auto.arima(resid, max.p=5, max.q=5)      # busqueda "automática" 
    cheq=Arima(resid, c(0,0,0), include.constant=TRUE) 
    summary(cheq)    
    tsdiag(cheq, gof=12)
  # Verificacion de normalidad en los residuales
    resid_sd=resid/ru_z@testreg$sigma
    qqnorm(resid_sd,  xlab = "Cuantiles Teóricos", ylab = "Cuantiles Muestrales")
    qqline(resid_sd)

    shapiro.test(resid_sd)             #prueba de Shapiro-Wilks
    jarqueberaTest(resid_sd )          # prueba Jarque-Bera 

# el modelo elegido para la prueba ADF satisface los supuestos básicos

# Conclusión: Parece que la serie z contiene al menos una raíz unitaria.

# Examen de más raíces unitarias
# gráfica de la serie diferenciada una vez
plot.ts(diff(z), type="o")
# prueba de si hay raíz unitarias en z diferenciada una vez.
ru_dif_z=ur.df(diff(z), type = c("drift"), lags=maxlag, selectlags = c("BIC"))
summary(ru_dif_z)  

# redefinición del modelo para la prueba
ru_dif_z=ur.df(diff(z), type = c("drift"), lags=0)
summary(ru_dif_z)  

# se validan y se elige el primer modelo
  # validación de la ecuación ADF
    resid1=ru_dif_z@testreg$residuals             # residuales del modelo ADF
    plot(ru_dif_z)                                # diagnósticos 
    auto.arima(resid1, max.p=5, max.q=5)       # busqueda "automática" 
    cheq1=Arima(resid1, c(0,0,0), include.constant=TRUE) 
    summary(cheq1)    
    tsdiag(cheq1, gof=12)
  # Verificacion de normalidad en los residuales
    qqnorm(resid1,  xlab = "Cuantiles Teóricos", ylab = "Cuantiles Muestrales")
    qqline(resid1)
    shapiro.test(resid1)                #prueba de Shapiro-Wilks
    library(fBasics)                   # librería fBasics
    jarqueberaTest(resid1)             # prueba Jarque-Bera 

# Conclusión: se cumplen los supuestos de la prueba. 
# La serie z diferenciada no contiene ráiz unitaria.
# Por tanto, la serie z pertenece a la clase de modelos ARIMA con deriva y d=1.

# determinación de los valores de (p, q) del modelo ARMA para (1-B)zt
# correlogramas muestrales para z diferenciada una vez
par(mfrow=c(2,1))
Acf(diff(z), lag.max=12, ci=0, ylim=c(-1,1))
Pacf(diff(z), lag.max=12, ci=0, ylim=c(-1,1))

# correlogramas muestrales  para (1-B)z con bandas
par(mfrow=c(2,1))
Acf(diff(z), lag.max=16)
Pacf(diff(z), lag.max=16)

# Parece que el proceso que generó la serie z es un ARIMA(0,1,0).
# selección usando criterios de información: use script criter_info_PIB_Col.r
# usando el criterio de información AIC indica el modelo ARIMA(0,1,3). 
# usando el criterio de información BIC indica el modelo ARIMA(0,1,0). 

# selección "automática" del modelo
auto.arima(z, d=1, max.p=5, max.q=5, ic=c("aic"))
auto.arima(z, d=1, max.p=5, max.q=5, ic=c("bic"))
# Usando el criterio de información AIC parece que el proceso que generó 
# la serie es un ARIMA(0,1,3) con deriva.
# Usando el criterio de información BIC parece que el proceso que generó 
# la serie es un ARIMA(0,1,0) con deriva.

# selección del modelo usando la eacf: se elige el modelo que señala el vértice 
# de un triángulo de ceros en la tabla                      
eacf(diff(z))
# Hay evidencia de que el proceso que generó la serie es un ARIMA(0,1,0) o un
# ARIMA(1,1,0).

# ETAPA DE ESTIMACIÓN
# ARIMA(0,1,0): estimación ML exacta del con valores iniciales dados por la estimación condicional
mod1_CSS_ML=Arima(z, c(0, 1, 0), include.drift=TRUE, lambda=1, method = c("CSS-ML"))
summary(mod1_CSS_ML) 
(res1_CSS_ML=residuals(mod1_CSS_ML))
res1_est=res1_CSS_ML/(mod1_CSS_ML$sigma2^.5)  # estandarización de los residuales

# ARIMA(1,1,0): estimación ML exacta del con valores iniciales dados por la estimación condicional
mod3_CSS_ML=Arima(z, c(1, 1, 0), include.drift=TRUE, lambda=1, method = c("CSS-ML"))
summary(mod3_CSS_ML) 
(res3_CSS_ML=residuals(mod3_CSS_ML))
res3_est=res3_CSS_ML/(mod3_CSS_ML$sigma2^.5)  # estandarización de los residuales

# ARIMA(0,1,3): estimación ML exacta del con valores iniciales dados por la estimación condicional
mod4_CSS_ML=Arima(z, c(0, 1, 3), include.drift=TRUE, lambda=1, method = c("CSS-ML"))
summary(mod4_CSS_ML) 
(res4_CSS_ML=residuals(mod4_CSS_ML))
res4_est=res4_CSS_ML/(mod4_CSS_ML$sigma2^.5)  # estandarización de los residuales

# ETAPA DE DIAGNÓSTICOS 
# raíces de los polinomios: 
# ARIMA(0,1,0) NO HAY POLINOMIOS AR ni MA.

# ARIMA(1,1,0)HAY POLINOMIO AR.
#raíces del polinomio AR
autoplot(mod3_CSS_ML)  

# ARIMA(0,1,3)HAY POLINOMIO MA.
#raíces del polinomio MA
autoplot(mod4_CSS_ML)  

# Conclusión: Parece que el ARIMA(1,1,0) no es adecuados.
# se continua con los modelos ARIMA(0,1,0) con deriva y ARIMA(0,1,3) con deriva.
# Por tanto la tendencia observada en la gráfica de z es una mezcla de
# tendencia aleatoria y deteminística lineal.

# Análisis de los residuales
# ARIMA(0,1,0)
checkresiduals(mod1_CSS_ML, lag=15)   # de la librería forecast
tsdiag(mod1_CSS_ML, gof=15)           # de la librería stats
# parece que el modelo deja dudas en explicar bien toda la estructura de dependencia
# aquí se detendrían los diagnósticos para ese modelo y se trataría de identificar un nuevo modelo.

# ARIMA(1,1,0)
tsdiag(mod3_CSS_ML, gof=15)

# ARIMA(0,1,3)
checkresiduals(mod4_CSS_ML, lag=15)   # de la librería forecast
tsdiag(mod4_CSS_ML, gof=15)
# No se rechaza que los residuales del modelo son ruido blanco 

# chequeo de normalidad
# gráfico cuantil-cuantil
qqnorm(res4_est,  xlab = "Cuantiles Teóricos", ylab = "Cuantiles Muestrales")
qqline(res4_est)
 
# histograma, densidad kernel y gráfico normal
plot(density(res4_est))
mu<-mean(res4_est)
sigm<-sd(res4_est)
x<-seq(-4,4,length=500)
y<-dnorm(x,mu,sigm)
hist(res4_est, prob=T, col="yellow")
lines(density(res4_est))
lines(x,y,lwd=2,col="blue")
# conclusión: No se detecta alejamiento fuerte de la normalidad

# pruebas de normalidad
shapiro.test(res4_est)                 # prueba de Shapiro-Wilks

# se usa la librería fBasics para la descripción de la distribución de los residuales, incluyendo
# los coeficientes de asimetría y curtosis muestrales y para la prueba de ljung_Box

basicStats(res4_est)  # kurtosis corresponde al "exceso de curtosis=curtosis-3" 
normalTest(res4_est, method=("jb"))  # librería fBasics: puede realizar otras pruebas
                                     # "ks" for the Kolmogorov-Smirnov one–sample test, 
                                     # "sw" for the Shapiro-Wilk test, 
                                     # "jb" for the Jarque-Bera Test, 
                                     # "da" for the D'Agostino Test. The default value is "ks"

# conclusión: No se rechaza la normalidad, a un nivel de significancia de alpha=0.05.
# Sin embargo, si se usa alpha>=0.09, la prueba de Shapiro-Wilk rechazaría normalidad. 

# detección de observaciones atípicas distantes
# chequeo de observaciones atípicas extremas (no es un análisis completo de outliers)
plot.ts(res4_est, type="o", ylim=c(-4,4))
abline(a=-3, b=0, col="red", lty=2)
abline(a=3, b=0, col="red", lty=2)

ind=(abs(res4_est)>3)
sum(ind)
(grupo=cbind(seq(z),res4_est, ind))

# valores ajustados del modelo para z transformada
(ajust=mod4_CSS_ML$fitted)     
# gráfico para los valores ajustados y los valores observados
ts.plot(z,ajust)   # gráfico de las series contra el tiempo
lines(z, col="black", type="o", cex=.5)
lines(ajust, col="red", type="o", cex=.5)

# gráfico de dispersión de valores observados vs valores ajustados
plot(as.vector(z),as.vector(ajust), type="p")   # gráfico de dispersión de la serie observada 
abline(0,1, col="red")                          # contra la serie ajustada

# Evaluación de la significancia estadística de los coeficientes estimados
# Construcción de los estadísticos de prueba, t
# estadísticos t
(t=coef(mod4_CSS_ML)/(diag(vcov(mod4_CSS_ML)))^.5)
# valores P
(val_p=2*pnorm(t, lower.tail=FALSE))
# a un nivel de significancia(aproximado) de alpha=0.05, el coficiente ma1 no es sigificativo;
# todos los demás sí lo son.

# Conclusión: el modelo pasa los diagnósticos.
# Es el modelo útil para pronosticar?

# =================================== EVALUACIÓN DE LOS PRONÓSTICOS

# se contruye el modelo usando 56 datos (conjunto de entrenamiento) y se dejan
# los últimos 8 (conjunto de prueba) para evaluar la capacidad de pronóstico del modelo 
# estimación ML exacta con valores iniciales dados por la estimación condicional
(mod_Evalpron=Arima(z[1:56], c(0, 1, 3), include.drift=T, method = c("CSS-ML")))
summary(mod_Evalpron) # compare con el modelo con todos los datos
res_eval=residuals(mod_Evalpron)
tsdiag(mod_Evalpron, gof=14)
checkresiduals(mod_Evalpron)
qqnorm(residuals(mod_Evalpron))
qqline(residuals(mod_Evalpron))
shapiro.test(residuals(mod_Evalpron))

(z_pred=forecast(mod_Evalpron, h=8, level=c(80, 95), fan=FALSE))

# gráfico de los pronósticos
plot(z_pred)

# o, alternativamente,
autoplot(z_pred)

# Evaluación de los prónosticos
# lista de los valores reales y los pronósticos
cbind(z[57:64], z_pred$mean)

ts.plot(z[57:64], z_pred$mean, z_pred$lower[,2], z_pred$upper[,2], type="o")
lines(z_pred$mean, col="red")
lines(z_pred$lower[,2], col="blue")
lines(z_pred$upper[,2], col="blue")

# precisión de los pronósticos
# --------------------------------------------------
(recm=(mean((z[57:64]-ts(z_pred$mean))^2))^.5)         # raíz cuadrada error cuadrático medio 
(recmp=100*(mean(((z[57:64]-ts(z_pred$mean))/z[57:64])^2))^.5) # raíz cuadrada error cuadrático medio en porcentaje

(eam=mean(abs(z[57:64]-ts(z_pred$mean))))         # error absoluto medio
(eamp=100*mean(abs((z[57:64]-ts(z_pred$mean))/z[57:64]))) # error absoluto medio en porcentaje

# descomposición del error cuadrático medio
# --------------------------------------------------
# error cuadrático medio
(ecm=mean((z[57:64]-ts(z_pred$mean))^2)) 

# proporción de sesgo en media
(prop_sesgom=(mean(z_pred$mean)-mean(z[57:64]))^2/ecm) 

# cálculo de varianza sesgadas
(num_pron=length(z_pred$mean))
(sigmap=(((num_pron-1)/num_pron)*var(z_pred$mean))^.5) # desv. estánd de los pronósticosdatos 
(sigmar=(((num_pron-1)/num_pron)*var(z[57:64]))^.5)  # desv. estánd de los datos

# proporción de sesgo en varianza
(prop_var=(sigmap-sigmar)^2/ecm)         

# proporción de covarianza
(prop_covar=2*(1-cor(z_pred$mean, z[57:64]))*sigmap*sigmar/ecm)

(prop_sesgom+prop_var+prop_covar)  # Chequeo

# EJERCICIO: Construya un modelo para el PIB usando la transformación logarítmica y compare
# sus resultados con el modelo anterior.

# Verdaderos Pronósticos
# --------------------------------------------------
par(mfrow=c(1,2))
(z_pron<-forecast(mod4_CSS_ML, h=8, level=c(80,95),  fan=F))
plot(z_pron)

(z_pron1<-forecast(mod4_CSS_ML, h=8, level=c(80,95),  fan=T))
plot(forecast(z_pron))

# o, alternativamente
autoplot(forecast(z_pron))
autoplot(forecast(z_pron1))

