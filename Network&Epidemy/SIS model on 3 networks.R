
################
# Main Program #
################

library(igraph)
library(ggplot2)
library(deSolve)
library(reshape2)
library(R.matlab)
library(Matrix)
library(poweRlaw)


#############
# Parametri #
#############

N = 100 #numero nodi
simlength <- 300  #giorni di simulazione epidemia
simlength2 <- 300
beta <- 0.5  #probabilità di trasmissione virus
mu <- 0.2  #probabilità che un nodo infetto diventi suscettibile


#----------------------------------------------------------------------------------------
#Modello SIS con eq differenziali#
#----------------------------------------------------------------------------------------

time=seq(from=1,to=simlength,by=1)
time2=seq(from=1,to=simlength2,by=1)

initial_state_values=c(S=N-1,I=1)  # S=suscettibili, I=infetti 
parameters=c(mu, beta)  #metto le probabilità in un unico vettore 'parameters'


sis_model <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    
    lambda=beta*(I/N) 
    dS=-lambda*S+mu*I
    dI=lambda*S-mu*I
    
    return(list(c(dS,dI)))
  }
  )
}

output<-as.data.frame(ode(y=initial_state_values, func=sis_model, parms=parameters, times=time))
out_long=melt(output,id="time")
ggplot(data = out_long,          
       aes(x = time, y = value/N, colour = variable, group = variable)) +  
  geom_line() +xlab("Tempo(giorni)")+ylab("Porzione della popolazione")+scale_color_discrete(name="Stato")





#----------------------------------------------------------------------------------------
#SCALE FREE#
#----------------------------------------------------------------------------------------

#funzione per generare networks di tipo Barabasi#

generate.network.B <- function(N,links.per.step){
  L <- matrix(nrow=0,ncol=2)
  deg <- integer(N)
  for (i in 2:N) {
    n.new <- min(links.per.step,i-1)
    linkto <- sample(i-1,n.new,prob=deg[1:(i-1)]+1)
    newlinks <- cbind(rep(i,n.new),linkto)
    L <- rbind(L,newlinks)
    deg[i] = deg[i] + n.new
    deg[linkto] = deg[linkto]+1
  }
  colnames(L) <- NULL
  L
}


plot.spread <- TRUE  #parametro per visualizzare l'andamento epidemico su grafico

links <- generate.network.B(N,2)  #la rete viene generata

infected <- logical(N)  #infected è un vettore logico di dimensione N con True se il nodo è Infetto
patientzero <- sample(N,1)  #viene randomicamente scelto un nodo tra gli N totali
infected[patientzero] <- TRUE  #a quel nodo randomicamente scelto viene associato il valore True (è il paziente zero)

vettore_num_infetti <- integer(simlength)  #vettore la cui iesima componente contiente il numero di infetti nell'iesimo giorno di pandemia
vettore_num_suscettibili <- integer(simlength)  #vettore la cui iesima componente contiente il numero di susciettibili nell'iesimo giorno di pandemia

  network.i <- graph.edgelist(links,directed=FALSE)
  fixlayout <- layout.kamada.kawai(network.i)
  node.colour <- rep("SkyBlue2",N)
  node.colour[patientzero] <- "red"
  summary(network.i)
  Max_degree = (max(degree(network.i)))
  hist(degree(network.i),xlim = c(0,25), breaks = seq(0,Max_degree,1),xaxp=c(0,25,5), xlab="gradi", ylab = "Frequenza", main = "Distribuzione gradi scale free")
  abline(h=seq(0,40,2), col="gray", lty="dotted" )
  fit1 <-fit_power_law(degree(network.i), 1)
  plot(network.i,layout=fixlayout, main="Tempo = 0", vertex.color=node.colour)
  
  a_sf<-matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:(length(links)/2)) {
      if (links[j,1]==i){
        a_sf[i, links[j,2]]<-1
        a_sf[links[j,2], i]<-1
      }
    }
  }
  
  writeMat('A_sf', x=as.matrix(a_sf))

for (i in 1:simlength){
  
  print(i)
  
  #Fase infezione#
  discordant.links <- which(xor(infected[links[,1]],infected[links[,2]]))  #trova gli indici dei nodi che sono collegati e in cui un nodo è infetto e l'altro suscettibile
  transmit <- rbinom(length(discordant.links),1,beta)  #determina randomicamente quali dei discordant.links trasmettono il virus
  #aggiornamento del vettore di infezione
  transmitter.links <- discordant.links[transmit==1]
  nodes.of.transmitter.links <- unique(as.vector(links[transmitter.links,1:2]))  #mette entrambi i nodi del transmitter.links in un singolo vettore; 'unique' filtra le ripetizioni
  infected[nodes.of.transmitter.links] <- TRUE  #vengono messi entrambi i nodi (ora infetti) su TRUE

  #Fase cura#
  infected.links <- which(infected) #trova gli indici dei nodi infetti
  curable <- rbinom(length(infected.links),1,mu) #determina randomicamente quali degli infetti si curano
  #aggiornamento del vettore di infezione
  cured.links <- infected.links[curable==1]
  infected[cured.links] <- FALSE #vengono messi i nodi curati (ora suscettibili) su FALSE
  
  vettore_num_infetti[i] <- sum(infected)
  vettore_num_suscettibili[i] <- N-vettore_num_infetti[i]
  print(paste0('infetti:', vettore_num_infetti[i]))
  print(paste0('suscettibili:', vettore_num_suscettibili[i]))
  
  #Stampa grafo con infetti e sucettibili#
  # if (plot.spread){
  #   node.colour[infected] <- "red"
  #   plot(network.i,layout=fixlayout, main=paste("Tempo =", i), vertex.color=node.colour)
  # }
}


#----------------------------------------------------------------------------------------
#networks alla fine dell'evoluzione#
node.colour[infected] <- "red"
plot(network.i,layout=fixlayout, main=paste("Momento attuale"), vertex.color=node.colour)


#----------------------------------------------------------------------------------------
#stampa andamento epidemico simulato#
plot(vettore_num_infetti/N, type = 'l', main = 'Simulazione epidemia', xlab = 'tempo(giorni)', ylab = '% popolazione')
lines(vettore_num_suscettibili/N, type = 'l', col='green')
time_n=seq(from=25,to=simlength,by=1)
vettore_num_infetti_n=vettore_num_infetti[25:simlength]
linear_model<-lm(vettore_num_infetti_n/N ~ time_n)
abline(h=as.numeric(unlist(linear_model[[1]][1])), col='red')
#----------------------------------------------------------------------------------------
#confronto simulazione-studio teorico#

out_conf<-melt(data.frame(output, vettore_num_infetti, vettore_num_suscettibili), id="time")
ggplot(data = out_conf,
       aes(x = time, y = value/N, colour = variable, group = variable)) + geom_line() +
       xlab("Tempo(giorni)")+ylab("Porzione della popolazione")+scale_color_discrete(name="Stato")
      


#----------------------------------------------------------------------------------------
#EQUAZIONE 7.7#
#----------------------------------------------------------------------------------------

#studiamo come varia la probabilità di diventare infetto nel tempo

#il numero di Nodi N potrebbe essere troppo elevato per questo studio quindi è stato ridotto a N2
#volendo si può far capitare il paziente zero negli N2 nodi e togliere i restanti ma si andrebbe ad
#alterare la struttura del grafo per cui ho preferito costruire un nuovo network per lo studio in questione

# links2 <- generate.network.B(N2,2)
# network2.i <- graph_from_edgelist(links2)
# 
# #troviamo la matrice di connettività
# a_sf<-matrix(0, N2, N2)
# for (i in 1:N2) {
#   for (j in 1:(length(links2)/2)) {
#     if (links2[j,1]==i){
#       a_sf[i, links2[j,2]]<-1
#       a_sf[links2[j,2], i]<-1
#     }
#   }
# }
# 
# writeMat('A_sf', x=as.matrix(a_sf))

#risolvo e grafico le eq differenziali
# sis_IBMF_model <- function(time,state,parameters){
#   with(as.list(c(state,parameters)),{
#     lambda=beta/mu
#     dp.i=-p.i+lambda*(1-p.i)*rowSums(a*p.i)
#     return(list(dp.i))
#   }
#   )
# }

# output_IBMF<-as.data.frame(ode(y=initial_state_values_IBMF,func = sis_IBMF_model,parms=parameters,times = time))
# out_long_IBMF=melt(output_IBMF,id="time")
# ggplot(data = out_long_IBMF,          
#        aes(x = time, y = value, colour = variable, group = variable)) +  
#   geom_line() +xlab("Tempo(giorni)")+ylab("Prob che iesimo nodo sia infetto(t)")+scale_color_discrete(name="Stato")


# #----------------------------------------------------------------------------------------
# #studio soglia#
# 
# #dalla sezione precendente possiamo ottenere lo jacobiano e determinarne l'autovalore di modulo max
# lambda=beta/mu
# Jac <- matrix(0, N2, N2)
# for (i in 1:N2) {
#   for (j in 1:N2) {
#     if (i==j)
#       Jac[i,j]=1+lambda*a[i,j]
#     else
#       Jac[i,j]=lambda*a[i,j]
#   }
# }
# 
# Aut <- eigen(Jac) #trovo autovalori e autovettori di Jac
# AutMax1 <- max(abs(Aut$values))  #trovo l'autovalore di max modulo
# Soglia1 <- 1/AutMax1  #calcolo la soglia
# 
# #Metodo teorico
# k_max <- max(degree(network2.i))
# k_mean <- mean(degree(network2.i))
# ksq_mean <- mean((degree(network2.i))^2)
# AutMax2 <- max(c(k_max, (ksq_mean/k_mean)))
# Soglia2 <- 1/AutMax2
# 
# #calcolo R0 per modello teorico eq diff
# Ro=beta*N-mu
# Sendemic=mu/beta



#----------------------------------------------------------------------------------------
#GRAFO OMOGENEO#
#----------------------------------------------------------------------------------------

network_hom.i<-erdos.renyi.game(N, 1/25)
links_hom <- as_edgelist(network_hom.i)

infected_hom <- logical(N)  #infected è un vettore logico di dimensione N con True se il nodo è Infetto
patientzero_hom <- sample(N,1)  #viene randomicamente scelto un nodo tra gli N totali
infected_hom[patientzero_hom] <- TRUE  #a quel nodo randomicamente scelto viene associato il valore True (è il paziente zero)

vettore_num_infetti_hom <- integer(simlength)  #vettore la cui iesima componente contiente il numero di infetti nell'iesimo giorno di pandemia
vettore_num_suscettibili_hom <- integer(simlength)  #vettore la cui iesima componente contiente il numero di susciettibili nell'iesimo giorno di pandemia

  
  fixlayout <- layout.kamada.kawai(network_hom.i)
  node.colour <- rep("SkyBlue2",N)
  node.colour[patientzero_hom] <- "red"
  summary(network_hom.i)
  Max_degree_hom = (max(degree(network_hom.i)))
  hist(degree(network_hom.i),xlim = c(0,16), breaks = seq(0,Max_degree_hom,1), xaxp=c(0,16,4), xlab="gradi", ylab = "Frequenza", main = "Distribuzione gradi Omogeneo")
  abline(h=seq(0,20,1), col="gray", lty="dotted" )
  plot(network_hom.i,layout=fixlayout, main="Tempo = 0", vertex.color=node.colour, main="Grafo omogeneo")

  a_hom<-matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:(length(links_hom)/2)) {
      if (links_hom[j,1]==i){
        a_hom[i, links_hom[j,2]]<-1
        a_hom[links_hom[j,2], i]<-1
      }
    }
  }
  
  writeMat('A_hom', x=as.matrix(a_hom))
  

for (i in 1:simlength){
  
  print(i)
  
  #Fase infezione#
  discordant.links <- which(xor(infected_hom[links_hom[,1]],infected_hom[links_hom[,2]]))  #trova gli indici dei nodi che sono collegati e in cui un nodo è infetto e l'altro suscettibile
  transmit <- rbinom(length(discordant.links),1,beta)  #determina randomicamente quali dei discordant.links trasmettono il virus
  #aggiornamento del vettore di infezione
  transmitter.links <- discordant.links[transmit==1]
  nodes.of.transmitter.links <- unique(as.vector(links_hom[transmitter.links,1:2]))  #mette entrambi i nodi del transmitter.links in un singolo vettore; 'unique' filtra le ripetizioni
  infected_hom[nodes.of.transmitter.links] <- TRUE  #vengono messi entrambi i nodi (ora infetti) su TRUE
  
  #Fase cura#
  infected.links <- which(infected_hom) #trova gli indici dei nodi infetti
  curable <- rbinom(length(infected.links),1,mu) #determina randomicamente quali degli infetti si curano
  #aggiornamento del vettore di infezione
  cured.links <- infected.links[curable==1]
  infected_hom[cured.links] <- FALSE #vengono messi i nodi curati (ora suscettibili) su FALSE
  
  vettore_num_infetti_hom[i] <- sum(infected_hom)
  vettore_num_suscettibili_hom[i] <- N-vettore_num_infetti_hom[i]
  print(paste0('infetti:', vettore_num_infetti_hom[i]))
  print(paste0('suscettibili:', vettore_num_suscettibili_hom[i]))
  
}


#----------------------------------------------------------------------------------------
#stampa andamento epidemico simulato#
plot(vettore_num_infetti_hom/N, type = 'l', main = 'Simulazione epidemia grafo hom', xlab = 'tempo(giorni)', ylab = '% popolazione')
lines(vettore_num_suscettibili_hom/N, type = 'l', col='green')
vettore_num_infetti_hom_n=vettore_num_infetti_hom[25:simlength]
linear_model_hom<-lm(vettore_num_infetti_hom_n/N ~ time_n)
abline(h=as.numeric(unlist(linear_model_hom[[1]][1])), col='red')

#networks alla fine dell'evoluzione#
node.colour[infected] <- "red"
plot(network_hom.i,layout=fixlayout, main=paste("fase finale"), vertex.color=node.colour)


#----------------------------------------------------------------------------------------
#confronto simulazione-studio teorico#

out_conf_hom<-melt(data.frame(output, vettore_num_infetti, vettore_num_suscettibili, vettore_num_infetti_hom, vettore_num_suscettibili_hom), id="time")
ggplot(data = out_conf_hom,
       aes(x = time, y = value/N, colour = variable, group = variable)) + geom_line() +
  xlab("Tempo(giorni)")+ylab("Porzione della popolazione")+scale_color_discrete(name="Stato")

#----------------------------------------------------------------------------------------
#matrice connettività hom#

# network3.i <- erdos.renyi.game(N2, 1/12)
# links3 <- as_edgelist(network3.i)
# 
# a_hom<-matrix(0, N2, N2)
# for (i in 1:N2) {
#   for (j in 1:(length(links3)/2)) {
#     if (links3[j,1]==i){
#       a_hom[i, links3[j,2]]<-1
#       a_hom[links3[j,2], i]<-1
#     }
#   }
# }
# 
# writeMat('A_hom', x=as.matrix(a_hom))


# #studio soglia#
# lambda=beta/mu
# Jac.hom <- matrix(0, N2, N2)
# for (i in 1:N2) {
#   for (j in 1:N2) {
#     if (i==j)
#       Jac.hom[i,j]=1+lambda*a.hom[i,j]
#     else
#       Jac.hom[i,j]=lambda*a.hom[i,j]
#   }
# }
# 
# Aut.hom <- eigen(Jac.hom) #trovo autovalori e autovettori di Jac
# AutMax1.hom <- max(abs(Aut.hom$values))  #trovo l'autovalore di max modulo
# Soglia1.hom <- 1/AutMax1.hom  #calcolo la soglia
# 
# #Metodo teorico
# k_max.hom <- max(degree(network3.i))
# k_mean.hom <- mean(degree(network3.i))
# ksq_mean.hom <- mean((degree(network3.i))^2)
# AutMax2.hom <- max(c(k_max.hom, (ksq_mean.hom/k_mean.hom)))
# Soglia2.hom <- 1/AutMax2.hom




#----------------------------------------------------------------------------------------
#GRAFO SMALL-WORLD#
#----------------------------------------------------------------------------------------

#generare networks small world#

sm <- sample_smallworld(1, N, 2, 0.05)
links_sm <- as_edgelist(sm)

fixlayout <- layout.kamada.kawai(sm)
node.colour <- rep("SkyBlue2",N)
Max_degree_sm = (max(degree(sm)))
hist(degree(sm),xlim = c(0,16), breaks = seq(0,Max_degree_sm,1), xaxp=c(0,16,4), xlab="gradi", ylab = "Frequenza", main = "Distribuzione gradi Small world")
abline(h=seq(0,80,5), col="gray", lty="dotted" )
plot(sm,layout=fixlayout, main="Tempo = 0", vertex.color=node.colour)

a_sm<-matrix(0, N, N)
for (i in 1:N) {
  for (j in 1:(length(links_sm)/2)) {
    if (links_sm[j,1]==i){
      a_sm[i, links_sm[j,2]]<-1
      a_sm[links_sm[j,2], i]<-1
    }
  }
}

writeMat('A_sm', x=as.matrix(a_sm))

infected_sm <- logical(N)  #infected è un vettore logico di dimensione N con True se il nodo è Infetto
patientzero_sm <- sample(N,1)  #viene randomicamente scelto un nodo tra gli N totali
infected_sm[patientzero_sm] <- TRUE  #a quel nodo randomicamente scelto viene associato il valore True (è il paziente zero)

vettore_num_infetti_sm <- integer(simlength)  #vettore la cui iesima componente contiente il numero di infetti nell'iesimo giorno di pandemia
vettore_num_suscettibili_sm <- integer(simlength)  #vettore la cui iesima componente contiente il numero di susciettibili nell'iesimo giorno di pandemia


for (i in 1:simlength){
  
  print(i)
  
  #Fase infezione#
  discordant.links <- which(xor(infected_sm[links_sm[,1]],infected_sm[links_sm[,2]]))  #trova gli indici dei nodi che sono collegati e in cui un nodo è infetto e l'altro suscettibile
  transmit <- rbinom(length(discordant.links),1,beta)  #determina randomicamente quali dei discordant.links trasmettono il virus
  #aggiornamento del vettore di infezione
  transmitter.links <- discordant.links[transmit==1]
  nodes.of.transmitter.links <- unique(as.vector(links_sm[transmitter.links,1:2]))  #mette entrambi i nodi del transmitter.links in un singolo vettore; 'unique' filtra le ripetizioni
  infected_sm[nodes.of.transmitter.links] <- TRUE  #vengono messi entrambi i nodi (ora infetti) su TRUE
  
  #Fase cura#
  infected.links <- which(infected_sm) #trova gli indici dei nodi infetti
  curable <- rbinom(length(infected.links),1,mu) #determina randomicamente quali degli infetti si curano
  #aggiornamento del vettore di infezione
  cured.links <- infected.links[curable==1]
  infected_sm[cured.links] <- FALSE #vengono messi i nodi curati (ora suscettibili) su FALSE
  
  vettore_num_infetti_sm[i] <- sum(infected_sm)
  vettore_num_suscettibili_sm[i] <- N-vettore_num_infetti_sm[i]
  print(paste0('infetti:', vettore_num_infetti_sm[i]))
  print(paste0('suscettibili:', vettore_num_suscettibili_sm[i]))
  
}


#stampa andamento epidemico simulato#
plot(vettore_num_infetti_sm/N, type = 'l', main = 'Simulazione epidemia grafo Small-World', xlab = 'tempo(giorni)', ylab = '% popolazione')
lines(vettore_num_suscettibili_sm/N, type = 'l', col='green')
vettore_num_infetti_sm_n=vettore_num_infetti_sm[25:simlength]
linear_model_sm<-lm(vettore_num_infetti_sm_n/N ~ time_n)
abline(h=as.numeric(unlist(linear_model_sm[[1]][1])), col='red')

#networks alla fine dell'evoluzione#
node.colour[infected] <- "red"
plot(sm,layout=fixlayout, main=paste("fase finale"), vertex.color=node.colour)



time2=seq(from=1,to=simlength2,by=1)

#confronto simulazione-studio teorico all epidemy#
out_conf_hom<-melt(data.frame(output, vettore_num_infetti, vettore_num_suscettibili, vettore_num_infetti_hom, vettore_num_suscettibili_hom, vettore_num_infetti_sm, vettore_num_suscettibili_sm), id="time")
ggplot(data = out_conf_hom,
       aes(x = time, y = value/N, colour = variable, group = variable)) + geom_line() +
  xlab("Tempo(giorni)")+ylab("Porzione della popolazione")+scale_color_discrete(name="Stato")


#confronto simulazione-studio teorico hom - deterministic#
out_conf_hom<-melt(data.frame(output, vettore_num_infetti_hom, vettore_num_suscettibili_hom), id="time")
ggplot(data = out_conf_hom,
       aes(x = time, y = value/N, colour = variable, group = variable)) + geom_line() +
  xlab("Tempo(giorni)")+ylab("Porzione della popolazione")+scale_color_discrete(name="Stato")



#confronto simulazione-studio teorico scale free - sm#
out_conf_hom<-melt(data.frame(output, vettore_num_infetti, vettore_num_suscettibili, vettore_num_infetti_sm, vettore_num_suscettibili_sm), id="time")
ggplot(data = out_conf_hom,
       aes(x = time, y = value/N, colour = variable, group = variable)) + geom_line() +
  xlab("Tempo(giorni)")+ylab("Porzione della popolazione")+scale_color_discrete(name="Stato")


#troviamo la matrice di connettività
# network4.i <- sample_smallworld(1, N2, 2, 0.05)
# links4 <- as_edgelist(network4.i)
# 
# a_sm<-matrix(0, N2, N2)
# for (i in 1:N2) {
#   for (j in 1:(length(links4)/2)) {
#     if (links4[j,1]==i){
#       a_sm[i, links4[j,2]]<-1
#       a_sm[links4[j,2], i]<-1
#     }
#   }
# }
# 
# writeMat('A_sm', x=as.matrix(a_sm))


# #studio soglia#
# lambda=beta/mu
# Jac.sm <- matrix(0, N2, N2)
# for (i in 1:N2) {
#   for (j in 1:N2) {
#     if (i==j)
#       Jac.sm[i,j]=1+lambda*a.sm[i,j]
#     else
#       Jac.sm[i,j]=lambda*a.sm[i,j]
#   }
# }
# 
# Aut.sm <- eigen(Jac.sm) #trovo autovalori e autovettori di Jac
# AutMax1.sm <- max(abs(Aut.sm$values))  #trovo l'autovalore di max modulo
# Soglia1.sm <- 1/AutMax1.sm  #calcolo la soglia
# 
# #Metodo teorico
# k_max.sm <- max(degree(network4.i))
# k_mean.sm <- mean(degree(network4.i))
# ksq_mean.sm <- mean((degree(network4.i))^2)
# AutMax2.sm <- max(c(k_max.sm, (ksq_mean.sm/k_mean.sm)))
# Soglia2.sm <- 1/AutMax2.sm


