#first let's work with the ordination, we need to assign value for each of the category if needed 
ps<-readRDS("data/ps_noDuplicates.RDS")
metada<-sample_data(ps)
metada<-as.data.frame(metada)

#removing the behavior questions 
metada<-metada[,-grep("Q",colnames(metada))]

#LEt's check something here
test<-metada[,c("food.allergies","Glutintol")]


#Adding column gluten free 
metada$gluten<-metada$food.allergies
levels(metada$gluten)[levels(metada$gluten)=="GlutenFree"]<-"Gluten.int"
levels(metada$gluten)[levels(metada$gluten)=="Dairy_free"]<-"None"
levels(metada$gluten)[levels(metada$gluten)== "GlutenANDdairy_free"]<-"Gluten.int"
levels(metada$gluten)[levels(metada$gluten)== "other"]<-"None"
metada$gluten[is.na(metada$gluten)]<-"None"
levels(metada$gluten)[levels(metada$gluten)==""]<-NA
metada$dairy<-metada$food.allergies
levels(metada$dairy)<-c(NA,"Dairy_int","Dairy_int","None","None")
metada$dairy[is.na(metada$dairy)]<-"None"
levels(metada$dairy)[levels(metada$dairy)==""]<-NA

#Removing kit_zip
metada<-metada[,-which(colnames(metada) =="kit.zip")]



#Consolidating child gastro status
metada$affecting_stool_constitancy<-metada$dietary.symptoms
levels(metada$affecting_stool_constitancy)<-c("","","","constipation","diarrhea","diarrhea/constipation","")
#now need to add the siblings
pair_nb_for_siblings_with_stool_cons<-metada$Pair[!is.na(metada$dietary.symptoms.sibling) & metada$Treatment == "Aut" ]
#check we have the corresponding siblings
pair_nb_for_siblings_with_stool_cons_cont<-metada$Pair[metada$Pair %in% pair_nb_for_siblings_with_stool_cons & metada$Treatment == "Control" ]
#looks good
metada$affecting_stool_constitancy[metada$Pair %in% pair_nb_for_siblings_with_stool_cons & metada$Treatment == "Control"]
metada$affecting_stool_constitancy<-as.character(metada$affecting_stool_constitancy)
metada$dietary.symptoms.sibling<-as.character(metada$dietary.symptoms.sibling)
metada$affecting_stool_constitancy[metada$Pair %in% pair_nb_for_siblings_with_stool_cons & metada$Treatment == "Control"]<-metada$dietary.symptoms.sibling[!is.na(metada$dietary.symptoms.sibling) & metada$Treatment == "Aut"]
metada$affecting_stool_constitancy<-as.factor(metada$affecting_stool_constitancy)
levels(metada$affecting_stool_constitancy)<-c("","constipation","diarrhea","","","")

#Now remove useless 4 columns 
metada<-metada[,-which(colnames(metada) %in% c("Child.Diagnosis","child.gastro","dietary.symptoms","pb.gastro.sibling","dietary.symptoms.sibling"))]

#one of the pair has "levels(metada$what.type)<-c("","Albendazole","Amoxicillin","Amoxicillin","Cefdinir")" => antiparasite! Albendazole
#some are several months ago, but onl 3 are this week 
#replacing with values: number of weeks: 
levels(metada$recent.antibiotics)<-c("",4,1,12,8)

#Reoving food allergy since we dealed with it earlier
metada<-metada[,-which(colnames(metada)  == "food.allergies")]

#also need to format this better 
levels(metada$ParentOtherInfo)<-c("Ab.pregnancy","gestational_diabetes",NA,"Preclampsia")
levels(metada$Probiotic)<-c(0,1,2,3,4)
levels(metada$VIT.B)<-c(0,1,2,3,4)
levels(metada$Vit.D)<-c(0,1,2,3,4)

#REmoving the lactose intolerent/Gluten intolerent questions since seems like there was a problem, see email to Jessica dn also it's reduncant with the other quetsions
metada<-metada[,-which(colnames(metada) %in% c("LactIntol","Glutintol"))]

#metada$WaterSou
#need to change the not sure into NA 
levels(metada$WaterSou)[6]<-"F"
levels(metada$WaterSou)[5]<-NA
levels(metada$teeth.Br)<-c(1,2,3,4)

#we have a "B/C" so we will put 1.5
levels(metada$bowel.freq)<-c(0,1,1.5,2,3,4)
levels(metada$Gender)<-c("f","f","m","m")
levels(metada$B.feed)[4]<-NA
levels(metada$Meat.eggs)<-c(0,1,2,3,4)
levels(metada$WholeGrain)<-c(0,1,2,3,4)
levels(metada$Fruit)<-c(0,1,2,3,4)
levels(metada$Veg)<-c(0,1,2,3,4)
levels(metada$Fermented.veg)<-c(0,1,2,3,4)
levels(metada$Milk..Cheese)<-c(0,1,2,3,4)
levels(metada$MilkSub)<-c(0,1,2,3,4)
levels(metada$Frozen.D)<-c(0,1,2,3,4)
levels(metada$Red.Meat)<-c(0,1,2,3,4)
levels(metada$Highfat.red.meat)<-c(0,1,2,3)
levels(metada$ChickenTurkey)<-c(0,1,2,3,4)
levels(metada$Seafood)<-c(0,1,2,3,4)
levels(metada$SaltedSnacks)<-c(0,1,2,3,4)
levels(metada$Sugar)<-c(0,1,2,3,4)
levels(metada$WholeEggs)<-c(0,1,2,3,4)
levels(metada$SweetBev)<-c(0,1,2,3,4)
levels(metada$Water.Day)<-c(0,1,2,3,4)
levels(metada$Ethnicity) #nothing to change here 
levels(metada$Current.Res)<-c(1,2,3,4,5)
levels(metada$Travel)<-c(NA,"C","D","E")
levels(metada$Cat)<-c("A","A","B","None")
levels(metada$OliveOil)<-c(0,1,2,3,4)
levels(metada$Dog)
levels(metada$bowel.qual)<-c("A","C","B","D")
#and then als change the age for the years
metada$age<-metada$age_month_ok
metada$age<-metada$age/12
metada$age<-ceiling(metada$age)
metada$age<-as.integer(metada$age)
metada<-metada[,-which(colnames(metada)=="age_month_ok")]
levels(metada$Gender)<-c("f","f","m","m")

#one last check
metada$Gender<-as.factor(metada$Gender)
colnames(metada)[which(colnames(metada) == "affecting_stool_constitancy")]<-"affecting_stool"
##Finally categorize the data from the metada$bowel.qual in order to have diahree and constip
metada$diarrhea<-metada$bowel.qual
levels(metada$diarrhea)<-c("other","constip","other","other")
metada$constip<-metada$bowel.qual
levels(metada$constip)<-c("diarrhea","other","other","other")
metada$gastro_pb<-metada$bowel.qual
levels(metada$gastro_pb)<-c("GI_pb","GI_pb","normal",NA)


save(metada,file="Metadata_microbiome_analysis/metada.rds")




