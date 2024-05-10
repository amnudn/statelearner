zel_sim2_1_n_events_mc_fun <- function(){
    ## Write sim setup here
    eval_times = seq(6, 36, 6)
    sim_sets = list(original = simZelefsky_wrapper,
                    indep_cens = simZelefsky_indep_cens_wrapper)
    out = do.call(rbind, lapply(seq_along(sim_sets), function(pp){
        sim_set_name = names(sim_sets)[pp]
        mc_sim = sim_sets[[pp]](n = 100000)
        do.call(rbind, lapply(eval_times, function(tt){
            data.table(sim_setting = sim_set_name,
                       time = tt,
                       true_events = mc_sim[, sum(true_time<tt)],
                       true_cens = mc_sim[, sum(cens_time<tt)],
                       at_risk = mc_sim[, sum(time>tt)])[]
        }))
    }))
    return(out[])
}
