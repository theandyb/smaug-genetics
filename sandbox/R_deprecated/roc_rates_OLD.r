require(ggplot2)
require(dplyr)
require(tidyr)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##############################################################################
# Read and process data
##############################################################################
# chrp<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/full/chr18_comb.txt", header=F)
# chrp<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/chr18_1pct_mask.txt", header=F)

# Read data
cat("Reading data...\n")
chrpf<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/full/rocdat_comb_7bp.txt", header=F)
names(chrpf)<-c("CHR", "POS", "MU", "OBS", "SEQ3", "MU3")

# Remove CpGs and sites with mu=0
chrpf<-chrpf[substr(chrpf$SEQ3, 2, 3)!="CG" & chrpf$MU>0,]

# Read DNMs
cat("Reading DNMs...\n")
dnms_full<-read.table("/net/bipolar/jedidiah/mutation/reference_data/DNMs/GoNL_DNMs.txt", header=T, stringsAsFactors=F)
dnms_full<-dnms_full[,1:3]
names(dnms_full)[1:3]<-c("ID", "CHR", "POS")

# Duplicate data, merge with DNMs to get ID
# chrpf<-chrp
cat("Annotating with ID...\n")
chrpf<-merge(chrpf, dnms_full, by=c("CHR", "POS"), all.x=T)
chrpf$ID[is.na(chrpf$ID)] <- "all"

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")
chrpfa<-chrpf[chrpf$ID=="all",]
chrpfa<-chrpfa[sample(nrow(chrpfa), 1000000),]
chrpfdnm<-chrpf[chrpf$ID!="all",]

# Combine data
cat("Creating combined data...\n")
chrp<-rbind(chrpfdnm, chrpfa) %>% arrange(MU)
chrp$prop <- cumsum(chrp$OBS)/sum(chrp$OBS)

##############################################################################
# Calculate individual-level AUC from the GoNL DNMs
#
# Each iteration calculates the AUC for the following:
# AUC under 3-mer model
# Permuted AUC under 3-mer model
# AUC under logit model
# Permuted AUC under logit model
##############################################################################
nperm<-258

aucperm<-rep(0,nperm)
aucperm3<-aucperm

ids<-unique(chrpfdnm$ID)
numind<-length(ids)
aucind<-rep(0, numind)
aucind3<-aucind

for(i in 1:nperm){
	### Run permutations
	cat("Permuting AUC", "(", i, "of", nperm, ")...\n")
	ndnms<-round(rnorm(1, 42.7, 10.3), 0)
	nsamp<-1000000

	chrpsub<-rbind(chrpfdnm[sample(nrow(chrpfdnm), ndnms),], chrpfa) %>%
		arrange(MU)
	chrpsub$prop <- cumsum(chrpsub$OBS)/sum(chrpsub$OBS)
	chrpsub2<-chrpsub[sample(nrow(chrpsub), nsamp),] %>%
	  arrange(MU) %>%
	  mutate(ntile=ntile(MU, 1000))

	chrpsub3<-chrpsub %>% arrange(MU3)
	chrpsub3$prop <- cumsum(chrpsub3$OBS)/sum(chrpsub3$OBS)
	chrpsub3a<-chrpsub3[sample(nrow(chrpsub3), nsamp),] %>%
		arrange(MU3, prop) %>%
		mutate(ntile=ntile(MU3, 1000))

	auctmp <- chrpsub2 %>% summarise(AUC=1-sum(prop)/nsamp)
	aucperm[i] <- auctmp$AUC
	auctmp3 <- chrpsub3a %>% summarise(AUC=1-sum(prop)/nsamp)
	aucperm3[i] <- auctmp3$AUC

	curid<-ids[i]
	### Get empirical AUC
	cat("Calculating AUC for", curid, "(", i, "of", numind, ")...\n")
	datid<-chrpfdnm %>% filter(ID==curid)

	tmpdat<-rbind(datid, chrpfa) %>% arrange(MU)
	tmpdat$prop<-cumsum(tmpdat$OBS)/sum(tmpdat$OBS)
	tmpdat2<-tmpdat[sample(nrow(tmpdat), nsamp),] %>%
	  arrange(MU) %>%
	  mutate(ntile=ntile(MU, 1000))

	### Repeat with 3-mer model
	tmpdat3<-tmpdat %>% arrange(MU3)
	tmpdat3$prop<-cumsum(tmpdat3$OBS)/sum(tmpdat3$OBS)
	tmpdat3s<-tmpdat3[sample(nrow(tmpdat3), nsamp),] %>%
		arrange(MU3, prop) %>%
		mutate(ntile=ntile(MU3, 1000))

	tmpauc <- tmpdat2 %>% summarise(AUC=1-sum(prop)/nsamp)
	aucind[i]<-tmpauc$AUC
	tmpauc3 <- tmpdat3s %>% summarise(AUC=1-sum(prop)/nsamp)
	aucind3[i]<-tmpauc3$AUC
}

##############################################################################
# Calculate AUC for each individual under the logit and 3-mer models
##############################################################################
# ids<-unique(chrpfdnm$ID)
# numind<-length(ids)
# aucind<-rep(0, numind)
# aucind3<-aucind
#
# for(i in 1:numind){
# 	curid<-ids[i]
# 	cat("Calculating AUC for", curid, "(", i, "of", numind, ")...\n")
# 	datid<-chrpfdnm %>% filter(ID==curid)
# 	tmpdat<-rbind(datid, chrpfa) %>% arrange(MU)
# 	tmpdat$prop<-cumsum(tmpdat$OBS)/sum(tmpdat$OBS)
#
# 	nsamp<-50000
#
# 	tmpdat2<-tmpdat[sample(nrow(tmpdat), nsamp),] %>%
# 	  arrange(MU) %>%
# 	  mutate(ntile=ntile(MU, 1000))
#
# 	tmpauc <- tmpdat2 %>% summarise(AUC=1-sum(prop)/nsamp)
# 	aucind[i]<-tmpauc$AUC
#
# 	### Repeat with 3-mer model
# 	tmpdat3<-tmpdat %>% arrange(MU3)
# 	tmpdat3$prop<-cumsum(tmpdat3$OBS)/sum(tmpdat3$OBS)
#
# 	tmpdat3s<-tmpdat3[sample(nrow(tmpdat3), nsamp),] %>%
# 		arrange(MU3, prop) %>%
# 		mutate(ntile=ntile(MU3, 1000))
#
# 	tmpauc3 <- tmpdat3s %>% summarise(AUC=1-sum(prop)/nsamp)
# 	aucind3[i]<-tmpauc3$AUC
# }

##############################################################################
# Function checks if elements in a exist in b
# Output is binary vector of length same as b
##############################################################################
toBin <- function(a,b){
	as.numeric(is.element(b,a))
}

##############################################################################
# Simulate mutations by randomly selecting sites and checking if
# Bernoulli(mu)==1; loop continues until we have simulated the same number
# of sites as in observed data
##############################################################################
nsim<-200

for(i in 1:nsim){
	cat("Running simulation", i, "of", nsim, "...\n")
	nsample <- 20000
	mutated <- c()
	nsites<-round(rnorm(1, 42.7, 10.3), 0)

	while(length(mutated) < nsites){
		rowind <- sample(nrow(chrp), nsample)
		row <- chrp[rowind,]
		mu <- row$MU

		batch <- sapply(mu, function(x) rbinom(1,1,x))
		mutated <- c(mutated, row$POS[which(as.logical(batch))])
	}

	# chrp$TMP <- toBin(mutated[1:nsites], chrp$POS)
	mus<-toBin(mutated, chrp$POS)
	nsimi<-sum(mus)
	chrp$last<-cumsum(mus)/nsimi

	colnames(chrp)[ncol(chrp)] <- paste0("simprop", i)
}

##############################################################################
# Subsample 100k sites, resort by mu, add percentile column, and coerce
# from wide to long format
#
# Can delete chrp after this step
##############################################################################
cat("Subsampling data...\n")
nsamp<-100000

chrp2<-chrp[sample(nrow(chrp), nsamp),] %>%
  arrange(MU) %>%
  mutate(ntile=ntile(MU, 1000)) %>%
  gather(group, val, .dots=grep("prop", names(chrp)))

auc <- chrp2 %>% group_by(group) %>% summarise(AUC=1-sum(val)/nsamp)

cat("Cleaning up memory...\n")
gc()

##############################################################################
# Summarize simulated data for ROC plots
##############################################################################
chrp2_nt_mean<-chrp2 %>%
  group_by(group, ntile) %>%
  summarise(val=mean(val))

# Get subsets of true and simulated data
chrp2_nt_sim<-chrp2_nt_mean[chrp2_nt_mean$group!="prop",]
chrp2_nt_obs<-chrp2_nt_mean[chrp2_nt_mean$group=="prop",]

# Calculate lower/upper range from simulations
# (can switch commented summarise() lines to use CI)
chrp2_lb<-chrp2_nt_sim %>%
  group_by(ntile) %>%
  # summarise(lb=t.test(val)$conf.int[1]) %>%
	summarise(lb=max(val)) %>%
  mutate(group="lb") %>%
  select(group, ntile, lb)

chrp2_ub<-chrp2_nt_sim %>%
  group_by(ntile) %>%
  # summarise(ub=t.test(val)$conf.int[2]) %>%
	summarise(ub=min(val)) %>%
  mutate(group="ub") %>%
  select(group, ntile, ub)

# Add bounds to data frame
chrp2_nt_obs$lb <- chrp2_lb$lb
chrp2_nt_obs$ub <- chrp2_ub$ub

##############################################################################
# Combine AUC data to generate ROC plot
##############################################################################
full_auc_dat<-rbind(data.frame(chrp2_nt_obs[,1:3]), data.frame(chrp3b))
names(full_auc_dat)[1]<-"Model"
full_auc_dat$Model<-ifelse(full_auc_dat$Model=="prop", "Logit", "3-mer")
full_auc_bounds<-chrp2_nt_obs[,c(2:5)]

ggplot()+
  geom_line(data=full_auc_dat,
		aes(x=(1000-ntile)/1000, y=1-val, group=Model, colour=Model), size=1.2)+
  geom_abline(intercept=0, slope=1)+
  geom_ribbon(data=full_auc_bounds,
		aes(x=(1000-ntile)/1000, y=1-val, ymin=1-ub, ymax=1-lb), alpha=0.2)+
	scale_colour_manual(values=cbbPalette[2:3])+
  coord_cartesian(xlim=c(0,1))+
	xlab("False Positive Rate")+
	ylab("True Positive Rate")+
  theme_bw()+
	theme(
		# legend.position="none",
		axis.title.x=element_text(size=20),
		axis.title.y=element_text(size=20))
  # scale_y_continuous(limits=c(0,1))+
  # scale_x_continuous(limits=c(0,1000), labels=c())
ggsave("/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png", height=4, width=7.25)

##############################################################################
# Get AUC data for distribution plots
##############################################################################
auc$model<-"simulated"
auc$obs<-"simulated"

# aucind_df<-data.frame(AUC=aucind, gp="ind")
aucsim<-auc[-1,2:4]
auclogitobs<-data.frame(AUC=aucind, model="logit", obs="observed")
auc3merobs<-data.frame(AUC=aucind3, model="3-mer", obs="observed")
auclogitperm<-data.frame(AUC=aucperm, model="logit", obs="permuted")
auc3merperm<-data.frame(AUC=aucperm3, model="3-mer", obs="permuted")

# All AUCs
plotaucfull<-rbind(aucsim, auclogitobs, auc3merobs, auclogitperm, auc3merperm)
plotaucfull$model<-as.factor(plotaucfull$model)
plotaucfull$model<-relevel(plotaucfull$model, "3-mer")

# No simulations
plotaucnosim<-rbind(auclogitobs, auc3merobs, auclogitperm, auc3merperm)
plotaucnosim$model<-as.factor(plotaucnosim$model)
plotaucnosim$model<-relevel(plotaucnosim$model, "3-mer")

# Only observed data (no simulation or permutations)
plotaucobsonly<-rbind(auclogitobs, auc3merobs)
plotaucobsonly$model<-as.factor(plotaucobsonly$model)
plotaucobsonly$model<-relevel(plotaucobsonly$model, "3-mer")

# Observed + Simulated (no permutations)
plotaucnoperm<-rbind(aucsim, auclogitobs, auc3merobs)
plotaucnoperm$model<-as.factor(plotaucnoperm$model)
plotaucnoperm$model<-relevel(plotaucnoperm$model, "3-mer")

plotaucmean<-plotaucobsonly %>%
	group_by(model, obs) %>%
	summarise(AUC=mean(AUC))

ggplot()+
	geom_density(data=plotaucfull, aes(x=AUC, colour=model, linetype=obs))+
	geom_density(data=aucsim, aes(x=AUC), colour="black")+
	scale_colour_manual(values=cbbPalette[c(2,3,1)])+
	# geom_vline(data=plotaucmean,
	# 	aes(xintercept=AUC, colour=model, linetype=obs))+
	theme_classic()+
	guides(linetype=F)+
	theme(legend.position="bottom")
ggsave("/net/bipolar/jedidiah/mutation/images/auc_hist.png", height=4, width=4)

# Add paternal age
ages<-read.table("/net/bipolar/jedidiah/mutation/reference_data/gonl_fam_age.txt", header=T)
names(ages)<-c("ID", "nDNM", "FatherAge", "MotherAge", "Coverage")
plotaucobsonly<-rbind(auclogitobs, auc3merobs)
plotaucobsonly$ID<-c(ids, ids)
plotaucobsonly<-merge(plotaucobsonly, ages, by="ID")

##############################################################################
# Additional data summaries for paternal age, etc.
##############################################################################
auc_age_cor<-plotaucobsonly %>%
	group_by(model) %>%
	summarise(cor=cor(AUC, FatherAge, method="spearman"),
		cor.p=cor.test(AUC, FatherAge, method="spearman")$p.value,
		cornum=cor(AUC, nDNM, method="spearman"),
		cornum.p=cor.test(AUC, nDNM, method="spearman")$p.value)

fao<-plotaucobsonly %>%
	filter(model=="logit") %>%
	filter(AUC>0.669 | AUC<0.597) %>%
	mutate(quart=ifelse(AUC>0.669, "yes", "no")) %>%
	dplyr::select(quart, FatherAge) %>%
	spread(quart, FatherAge)

# OLD VERSION--early attempt at pseudo-ROC curves
# ggplot(chrp3, aes(x=1000-ntile, y=1-val, group=group, colour=group))+
#   geom_line()
# ggsave("/net/bipolar/jedidiah/mutation/images/psuedo_roc_chr4a.png")

# OLD VERSION--uses true ROC calculations
# devtools::install_github("hadley/ggplot2")
# devtools::install_github("sachsmc/plotROC")
# library(plotROC)
#
# basicplot2<-ggplot(chrp2, aes(d = OBS, m = MU))+
# 	geom_roc()
#
# basicplot2+annotate("text", x = .75, y = .25,
#            label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2)))+
# 	style_roc()
# ggsave("/net/bipolar/jedidiah/mutation/images/chr18_roc2.png")

############################################
# Server dump below here
############################################
# dnm_agg3 <- chrpfdnm %>%
# 	group_by(Category.x=as.character(Category.x), SEQ=substr(SEQ, 3, 5)) %>%
# 	summarise(ndnm=n())
#
# rates3 <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_1000k_rates.txt", header=T, stringsAsFactors=F)
#
# rates3$Category.x <- gsub("cpg_", "", rates3$Category2)
# rates3$SEQ <- substr(rates3$Sequence, 1, 3)
# r3m <- merge(dnm_agg, rates3, by=c("Category.x", "SEQ"))
#
# rates_3 <- r3m %>%
# 	mutate(rel_prop3=ndnm/COUNT) %>%
#   dplyr::select(Sequence, Category.x, rel_prop3) %>%
#   spread(Category.x, rel_prop3)
#
# rates_3[is.na(rates_3)] <- 0

chrp_maxc_d <- chrp_maxc %>%
  filter(!grepl("^goldmann", ID)) %>%
  group_by(Category.x, SEQ) %>%
  summarise(ndnms=sum(OBS)) %>%
  mutate(SEQ3=substr(SEQ, 3, 5),
    SEQ5=substr(SEQ, 2, 6),
    SEQ7=substr(SEQ, 1, 7))

getDNMRates <- function(mlength){
  infile <- paste0("/net/bipolar/jedidiah/mutation/output/", mlength, "bp_1000k_rates.txt")
  rates <- read.table(infile, header=T, stringsAsFactors=F)

  rates$Category.x <- gsub("cpg_", "", rates$Category2)
  rates$SEQ <- substr(rates$Sequence, 1, mlength)

  dnm_agg <- chrpfdnm %>%
    # filter(grepl("^goldmann", ID)) %>%
  	group_by(Category.x=as.character(Category.x),
    SEQ=substr(SEQ, 4-(mlength-1)/2, 4+(mlength-1)/2)) %>%
  	summarise(ndnm=n())

  r3m <- merge(dnm_agg, rates, by=c("Category.x", "SEQ"))

  rates_out <- r3m %>%
  	mutate(dnm_rel_prop=ndnm/COUNT,
      Category2=ifelse(substr(Sequence,4,5)=="CG",
        paste0("cpg_", Category.x),
        Category.x)) %>%
    filter(ndnm>=4) %>%
    dplyr::select(Sequence, SEQ, Category2, ndnm, dnm_rel_prop) #%>%
    # spread(Category.x, rel_prop)
}

r3a <- getDNMRates(3)
r5a <- getDNMRates(5)
r7a <- getDNMRates(7)

r2<-rates %>%
  group_by(Category2) %>%
  mutate(rk=ntile(desc(rel_prop), 100)) %>%
  filter(rk<21)
rates_m <- merge(rates, r7a, by=c("Category2", "Sequence"))
rmc <- rates_m %>%
  # filter(ndnm>=6) %>%
  group_by(Category2) %>%
  summarise(cor=cor(rel_prop, dnm_rel_prop),
    cor_s=cor(rel_prop, dnm_rel_prop, method="spearman"))

ggplot(rates_m, aes(x=dnm_rel_prop, y=rel_prop, colour=Category2))+
  geom_point()+
  facet_wrap(~Category2, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/dnm_rate_scatter.png")

chrp_maxc_d2 <- dplyr::select(
  merge(chrp_maxc_d, r3a,
    by.x=c("Category.x", "SEQ3"), by.y=c("Category.x", "SEQ")),
  Category.x, SEQ, SEQ3, SEQ5, SEQ7, ndnms, rp3=rel_prop)

chrp_maxc_d2 <- dplyr::select(
  merge(chrp_maxc_d2, r5a,
    by.x=c("Category.x", "SEQ5"), by.y=c("Category.x", "SEQ"), all.x=T),
  Category.x, SEQ, SEQ3, SEQ5, SEQ7, ndnms, rp3, rp5=rel_prop)

chrp_maxc_d2 <- dplyr::select(
  merge(chrp_maxc_d2, r7a,
    by.x=c("Category.x", "SEQ7"), by.y=c("Category.x", "SEQ"), all.x=T),
  Category.x, SEQ, SEQ3, SEQ5, SEQ7, ndnms, rp3, rp5, rp7=rel_prop)

chrp_maxc_d3 <- chrp_maxc_d2 %>%
  mutate(rel_prop=ifelse(ndnms>0, rp7, ifelse(!is.na(rp5), rp5, rp3))) %>%
  dplyr::select(Category.x, Sequence=SEQ, rel_prop) #%>%
  # spread(Category.x, rel_prop)

chrp_maxc_d3[is.na(chrp_maxc_d3)] <- 0


chrp_maxc_m <- merge(chrp_maxc[!grepl("^gonl", chrp_maxc$ID),], chrp_maxc_d3,
  by.x=c("Category.x", "SEQ"), by.y=c("Category.x", "Sequence"))

maxc_auc_dat2 <- subWrapper(chrp_maxc_m, sim=F)
chrp_m <- subByModel(chrp_maxc_m, "rel_prop", "DNMs")

maxc_auc_dat3 <- bind_rows(maxc_auc_dat2, chrp_m)

fig_max_dat2 <- maxc_auc_dat3 %>%
	filter(grepl("^E|^B|^D", group)) %>%
  filter(!grepl("^gonl", ID))
plotROC(fig_max_dat2, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves_maxc_dnm_2.png")
