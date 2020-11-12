setwd("C:/Users/lauren erdman/Desktop/pyloplasty/")

library(randomForest)
library(readxl)

pyl = data.frame(read_xlsx("Pyeloplasty.xlsx"))
head(pyl)


pyl$Anomalies[is.na(pyl$Anomalies)] = 0
pyl$Anomalies = factor(pyl$Anomalies,levels = 0:5)
pyl$Sex = factor(pyl$Sex, levels = c(1,2), labels = c("M","F"))
pyl$Surgeon = factor(pyl$Surgeon,levels = 1:9)
pyl$Side = factor(pyl$Side, levels = 1:2)
pyl$Who_indicated = factor(pyl$Who_indicated,1:5)
pyl$Approach = factor(pyl$Approach,1:4)
pyl$Intraop_finding[is.na(pyl$Intraop_finding)] = 0
pyl$Intraop_finding = factor(pyl$Intraop_finding, levels = 0:2)
pyl$Angiopexy = factor(pyl$Angiopexy, levels = 0:1)
pyl$Salle_Stent = factor(pyl$Salle_Stent, levels = 0:1)
pyl$JJ_Stent = factor(pyl$JJ_Stent, levels = 0:1)
pyl$Nx_Vx_Prophy = factor(pyl$Nx_Vx_Prophy, levels = 0:1)
pyl$Dex = factor(pyl$Dex, levels = 0:1)
pyl$Blocks = factor(pyl$Blocks, levels = 0:1)
pyl$NG_tube = factor(pyl$NG_tube, levels = 0:1)
pyl$Catheter = factor(pyl$Catheter, levels = 0:1)
pyl$Drain = factor(pyl$Drain, levels = 0:1)
pyl$Narc_Rx = factor(pyl$Narc_Rx, levels = 0:1)
pyl$PCA_use = factor(pyl$PCA_use, levels = 0:1)
pyl$Epidural = factor(pyl$Epidural, levels = 0:1)
pyl$ketoralac = factor(pyl$ketoralac, levels = 0:1)
pyl$IV_Fluid = factor(pyl$IV_Fluid, levels = 0:1)
pyl$Duration_IV[is.na(pyl$Duration_IV)] = 0
pyl$Constipation_prophy = factor(pyl$Constipation_prophy, levels = 0:1)
pyl$Mobilization = factor(pyl$Mobilization, levels = 0:1)
pyl$Return_Diet = factor(pyl$Return_Diet, levels = 0:4)

### time-variant train/val/test split

### OUTCOMES: emergency room (ER) visits/readmissions, 
###     unplanned additional procedures, and redo pyeloplasty
outcome_cols = c("Reoperation", "time_to_event_allmo", "REDO" ,"Bounceback" ,  "BB_Reason", 
                "Re_admit"  , "Preop_US",  "Preop_nuc_scan"   , "VCUG"  , 
                "RPG" ,  "ABdo_Xray"  , "Nephrostogram" ,"LOS_days" , "LOS_hours" )

### PREDICTORS: AGE, SEX, ABD (BASELINE AND FU), DRF
###       LOS, use of drains, catheters and/or stents, opioids, 
###       regional anesthesia, non-opioid analgesia, 
###       nausea/vomiting prophylaxis, return to diet, and mobilization

pred_cols = c("Year_Sx","Sex","Age_sx_Mos","Anomalies","Who_indicated"  ,    
              "Sex_Provider"      ,  "Surgeon"            ,
              "Side"              ,  "Approach"           , "Intraop_finding",
              "Angiopexy"         ,  "Salle_Stent"        , "JJ_Stent"  ,
              "OR_Time"           ,  "Nx_Vx_Prophy"       , "Dex" ,
              "Blocks"            ,
              "NG_tube"           ,   "Catheter"  ,
               "Drain"              , 
              "Narc_Rx"           ,  "Narc_use"           , 
              "PCA_use"           ,  "Epidural"  ,
               "ketoralac"          , "IV_Fluid"  ,
              "Duration_IV"       ,  "Constipation_prophy", "Mobilization"  ,
              "Return_Diet"       )


saveRDS(pyl, file = "C:/Users/lauren erdman/Desktop/pyloplasty/pyloplasty_preproc.rds")
