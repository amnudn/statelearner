### run.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2023 (14:07) 
## Version: 
## Last-Updated: Jul 23 2023 (19:48) 
##           By: Anders Munch
##     Update #: 173
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## library(here)
try(setwd("~/research/SuperVision/Anders/survival-loss/statelearner/empirical-study/"))
try(setwd("/home/amnudn/Documents/phd/survival-loss-function/statelearner/empirical-study"))
library(targets)
library(lava)
library(data.table)
library(riskRegression)
library(survival)
library(ranger)
library(ggplot2)
tar_make()

## IPCW fail example:
tar_load(ipcw_fail)
ipcw_fail[loss_type == "sl-brier", learner := "statelearner"]
ipcw_fail[loss_type == "brier", learner := paste0("IPCW(", cens_model , ")")]
winners <- rbind(ipcw_fail[!is.na(loss), .SD[min(loss) == loss, .(winner = out_model[1], oracle_loss = oracle_loss[1])], .(learner, n, n_sim)],
                 ipcw_fail[learner == "IPCW(dgm)", .SD[min(oracle_loss) == oracle_loss, .(winner = out_model[1], oracle_loss = oracle_loss[1])], .(learner, n, n_sim)][, learner := "oracle"])
winners_summ <- winners[, .(correct = mean(winner == "dgm")), .(learner,n)]

## Compare oracle risk
ggplot(winners[,.(oracle_risk = 100*mean(oracle_loss)), .(learner, n)], aes(x = n, y = oracle_risk, col = learner)) +
    geom_line() + geom_point() + 
    theme_bw()

## See how many times we get things right
ggplot(winners_summ, aes(x = n, y = correct, col = learner)) +
    geom_line() + geom_point() + 
    theme_bw() + 
    ylim(c(0,1))


## Simulate data based on zelefsky data
tar_load(sim_zel_learners)
z_out_winner <- do.call(rbind, lapply(sim_zel_learners[, unique(n_sim)], function(ii){
    w0 = rbind(sim_zel_learners[n_sim == ii & out_model != "pre-KM" & learner != "oracle",
                                .(winner = .SD[min(loss) == loss, out_model[1]]), .(learner, n)],
               sim_zel_learners[n_sim == ii & out_model != "pre-KM" & learner == "oracle" & cens_model == "no-cens",
                                .(winner = .SD[min(loss) == loss, out_model[1]]), .(learner, n)])
    oracle0 = sim_zel_learners[n_sim == ii & learner == "oracle" & cens_model == "no-cens", .(n, model = out_model, oracle_loss = loss)]
    out = merge(w0, oracle0, by.x = c("winner", "n"), by.y = c("model", "n"), all.x = TRUE)
    out[, oracle_loss := 100*oracle_loss]
    out[, n_sim := ii]
    return(out[])
}))
z_out_winner_summ <- z_out_winner[, .(loss = mean(oracle_loss, na.rm = TRUE), sd = sd(oracle_loss, na.rm = TRUE)), .(n,learner)]

## Oracle risk for outcome model
ggplot(z_out_winner_summ, aes(x = n, y = loss)) +
    geom_ribbon(aes(ymin = loss-2*sd/sqrt(200), ymax = loss+2*sd/sqrt(200), fill = learner), alpha = .3) +
    geom_line(aes(col = learner)) + geom_point(aes(col = learner)) +
    theme_bw()

z_out_winner[learner != "ipcw-km", .(agree = 1*(.SD[1, winner] == .SD[2, winner])), .(n,n_sim)][, mean(agree, na.rm = TRUE), n]
z_out_winner[learner != "sl", .(agree = 1*(.SD[1, winner] == .SD[2, winner])), .(n,n_sim)][, mean(agree, na.rm = TRUE), n]

z_cens_winner <- do.call(rbind, lapply(sim_zel_learners[, unique(n_sim)], function(ii){
    w0 = rbind(sim_zel_learners[n_sim == ii & cens_model != "pre-KM" & learner != "oracle",
                                .(winner = .SD[min(loss) == loss, cens_model[1]]), .(learner, n)],
               sim_zel_learners[n_sim == ii & cens_model != "pre-KM" & learner == "oracle" & out_model == "no-cens",
                          .(winner = .SD[min(loss) == loss, cens_model[1]]), .(learner, n)])
    oracle0 = sim_zel_learners[n_sim == ii & learner == "oracle" & out_model == "no-cens", .(n, model = cens_model, oracle_loss = loss)]
    out = merge(w0, oracle0, by.x = c("winner", "n"), by.y = c("model", "n"), all.x = TRUE)
    out[, oracle_loss := 100*oracle_loss]
    out[, n_sim := ii]
    return(out[])
}))
z_cens_winner_summ <- z_cens_winner[, .(loss = mean(oracle_loss, na.rm = TRUE), sd = sd(oracle_loss, na.rm = TRUE)), .(n,learner)]

## Oracle risk for outcome model
ggplot(z_cens_winner_summ, aes(x = n, y = loss)) +
    ## geom_ribbon(aes(ymin = loss-2*sd/sqrt(200), ymax = loss+2*sd/sqrt(200), fill = learner), alpha = .3) +
    geom_line(aes(col = learner)) + geom_point(aes(col = learner)) +
    theme_bw()





## ----------------- sandbox -----------------
z_cens_winner[learner != "ipcw-km", .(agree = 1*(.SD[1, winner] == .SD[2, winner])), .(n,n_sim)][, mean(agree, na.rm = TRUE), n]
z_cens_winner[learner != "sl", .(agree = 1*(.SD[1, winner] == .SD[2, winner])), .(n,n_sim)][, mean(agree, na.rm = TRUE), n]

sz_learners_out <- copy(sim_zel_learners)[out_model != "no-cens" & out_model != "pre-KM"]

sz_learners_out[learner == "oracle", loss_category := "oracle_loss"]
sz_learners_out[learner != "oracle", loss_category := "emp_loss"]
dcast(sz_learners_out, out_model + cens_model + learner + n + n_sim ~ loss_category, value.var = "loss")


z_winners_out <- sim_zel_learners[!is.na(loss) & out_model != "pre-KM", {
    w0 = do.call(rbind, lapply(.SD[out_model != "oracle", out_model], function(mm){        
        winner = .SD[min(loss) == loss, out_model[1]]
        winner
    }))
    w0
}, .(n, n_sim)]

merge(x = z_winners_out, y = z_oracle_loss_out,
      by.x = c("n_sim", "winner"), by.y = c("n_sim", "out_model"))



sim_zel_learners[learner != "oracle", .(n_sim, oracle_loss = loss)]

tar_load(zelefsky_statelearner)




tar_load(zelefsky_statelearner_real_data)
zelefsky_statelearner_real_data


library(ggplot2)
ggplot(zelefsky_statelearner, aes(x = out_model, y = loss)) +
    geom_errorbar(aes(ymin = loss-2*sd, ymax = loss+2*sd)) +
    geom_point() +
    facet_wrap(~cens_model) + theme_bw()

tar_load(test_fit)
test_fit
tar_load(test_ipcw)
test_ipcw

tar_load(simulated_zelefsky_data)
tar_load(zelefsky_fit)
Publish::publish(zelefsky_fit)

tar_load(simulation_input)

coxph(Surv(time,event)~logPSA*sDose + stage + ggtot+hormones, data=sim_z, x = TRUE, y = TRUE)

coxph(Surv(time,event)~pspline(logPSA) + sDose + stage + ggtot+hormones, data=sim_z, x = TRUE, y = TRUE)
coxph(Surv(time,event)~logPSA + pspline(sDose) + stage + ggtot+hormones, data=sim_z, x = TRUE, y = TRUE)
coxph(Surv(time,event)~logPSA + sDose + stage + pspline(ggtot)+hormones, data=sim_z, x = TRUE, y = TRUE)
coxph(Surv(time,event)~pspline(logPSA) + pspline(sDose) + stage + pspline(ggtot)+hormones, data=sim_z, x = TRUE, y = TRUE)

coxph(Surv(time,event == 0)~pspline(logPSA) + pspline(sDose) + stage + pspline(ggtot) + hormones, data=sim_z, x = TRUE, y = TRUE)

coxph(Surv(time,event == 0)~logPSA + pspline(sDose) + stage + ggtot + hormones, data=sim_z, x = TRUE, y = TRUE)

coxph(Surv(time,event == 0)~sDose + stage , data=sim_z, x = TRUE, y = TRUE)

sim_z <- simZelefsky(n = 2000,censoring = TRUE,simulation_input = simulation_input)
sim_z[, ":="(time = dmos, event = status)]
out_learners <- list(km = coxph(Surv(time,event)~1, data=sim_z, x = TRUE, y = TRUE),
                     cox_full = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=sim_z, x = TRUE, y = TRUE),
                     cox_int = coxph(Surv(time,event)~pspline(logPSA) + sDose + stage + ggtot+hormones, data=sim_z, x = TRUE, y = TRUE),
                     cox_spline = coxph(Surv(time,event)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_z, x = TRUE, y = TRUE),
                     cox_small = coxph(Surv(time,event)~logPSA+sDose+stage, data=sim_z, x = TRUE, y = TRUE)
                     )
x <- Score(out_learners, data = sim_z,formula = Hist(time,event)~1,metrics = "brier",se.fit = FALSE,times = 36,contrasts = FALSE,null.model = FALSE, B = 1, split.method = "cv5")
setorder(x$Brier$score, Brier)
x$Brier$score

mc_data <- simZelefsky(n = 20000,censoring = TRUE,simulation_input = simulation_input)
mc_data[, dummy := 1]
Score(out_learners, data = mc_data,formula = Hist(time.recur,dummy)~1,metrics = "brier",se.fit = FALSE,times = 36,contrasts = FALSE,null.model = FALSE)

Score(out_learners, data = sim_z,formula = Hist(time.recur,dummy)~1,metrics = "brier",se.fit = FALSE,times = 36,contrasts = FALSE,null.model = FALSE, B = 1, split.method = "cv5")$Brier$score

cens_learners <- list(km = coxph(Surv(time,event==0)~1, data=sim_z, x = TRUE, y = TRUE),
                      cox_full = coxph(Surv(time,event==0)~logPSA+stage+ggtot+sDose+hormones, data=sim_z, x = TRUE, y = TRUE),
                      cox_int = coxph(Surv(time,event==0)~pspline(logPSA) + sDose + stage + ggtot+hormones, data=sim_z, x = TRUE, y = TRUE),
                      cox_spline = coxph(Surv(time,event==0)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_z, x = TRUE, y = TRUE),
                      cox_small = coxph(Surv(time,event==0)~logPSA+sDose+stage, data=sim_z, x = TRUE, y = TRUE)
                      )

sim_z[, dummy := 1]
x_cens <- Score(cens_learners, data = sim_z,formula =Hist(time, event == 0)~logPSA+stage+ggtot+sDose+hormones ,metrics = "brier",se.fit = FALSE,times = 36,contrasts = FALSE,null.model = FALSE, B = 1, split.method = "cv5")
setorder(x_cens$Brier$score, Brier)
x_cens

x2 <- statelearner(learners = list(state = out_learners,
                                   censoring = cens_learners),
                   data = sim_z,
                   B = 1,
                   times = 36)
x2


######################################################################
### run.R ends here
