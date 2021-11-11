library(ape)

#test prefixes
#test_prefix_4a <- "test_4a"
#test_prefix_amb <- "test_amb"
#test_prefix_n <- "test_n"
#test_prefix_no_change <- "test_no_change"
#test_prefix_lb_1 <- "test_lb_1"
#test_prefix_lb_5 <- "test_lb_5"
#test_prefix_lb_10 <- "test_lb_10"

#tree_list <- list("test_4a", "test_amb", "test_n", "test_no_change", "test_lb_1", "test_lb_5", "test_lb_10")
#for(test_prefix in (tree_list)) {
 #   input_tree <- read.tree("/project/exaptation/fast_approx_pruning/test_data/") + test_prefix + "/" + ".fastree"
  #  bifurcating_tree <- multi2di(input_tree)
   # write.tree(bifurcating_tree,file = "/project/exaptation/fast_approx_pruning/test_data/") + test_prefix + "/" + "_bi.fastree"
    #}
    
    
    
    
input_tree_4a <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_4a/test_4a.fastree")
input_tree_n <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_n/test_n.fastree")
input_tree_amb <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_amb/test_amb.fastree")
input_tree_lb_1_perc <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_lb_1_perc/test_lb_1_perc.fastree")
input_tree_lb_5_perc <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_lb_5_perc/test_lb_5_perc.fastree")
input_tree_lb_10_perc <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_lb_10_perc/test_lb_10_perc.fastree")
input_tree_lb_1_change <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_lb_1_change/test_lb_1_change.fastree")
input_tree_lb_5_change <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_lb_5_change/test_lb_5_change.fastree")
input_tree_lb_10_change <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_lb_10_change/test_lb_10_change.fastree")
input_tree_no_change <- read.tree("/project/exaptation/fast_approx_pruning/test_data/test_no_change/test_no_change.fastree")

#is.rooted(input_tree)

# root the unrooted tree
bifurcating_tree_4a <- multi2di(input_tree_4a)
bifurcating_tree_n <- multi2di(input_tree_n)
bifurcating_tree_amb <- multi2di(input_tree_amb)
bifurcating_tree_lb_1_perc <- multi2di(input_tree_lb_1_perc)
bifurcating_tree_lb_5_perc <- multi2di(input_tree_lb_5_perc)
bifurcating_tree_lb_10_perc <- multi2di(input_tree_lb_10_perc)
bifurcating_tree_lb_1_change <- multi2di(input_tree_lb_1_change)
bifurcating_tree_lb_5_change <- multi2di(input_tree_lb_5_change)
bifurcating_tree_lb_10_change <- multi2di(input_tree_lb_10_change)
bifurcating_tree_no_change <- multi2di(input_tree_no_change)

# write tree to file
write.tree(bifurcating_tree_4a,file = "/project/exaptation/fast_approx_pruning/test_data/test_4a/test_4a_bi.fastree")
write.tree(bifurcating_tree_n,file = "/project/exaptation/fast_approx_pruning/test_data/test_n/test_n_bi.fastree")
write.tree(bifurcating_tree_amb,file = "/project/exaptation/fast_approx_pruning/test_data/test_amb/test_amb_bi.fastree")
write.tree(bifurcating_tree_lb_1_perc,file = "/project/exaptation/fast_approx_pruning/test_data/test_lb_1_perc/test_lb_1_perc_bi.fastree")
write.tree(bifurcating_tree_lb_5_perc,file = "/project/exaptation/fast_approx_pruning/test_data/test_lb_5_perc/test_lb_5_perc_bi.fastree")
write.tree(bifurcating_tree_lb_10_perc,file = "/project/exaptation/fast_approx_pruning/test_data/test_lb_10_perc/test_lb_10_perc_bi.fastree")
write.tree(bifurcating_tree_no_change,file = "/project/exaptation/fast_approx_pruning/test_data/test_no_change/test_no_change_bi.fastree")
write.tree(bifurcating_tree_lb_1_change,file = "/project/exaptation/fast_approx_pruning/test_data/test_lb_1_change/test_lb_1_change_bi.fastree")
write.tree(bifurcating_tree_lb_5_change,file = "/project/exaptation/fast_approx_pruning/test_data/test_lb_5_change/test_lb_5_change_bi.fastree")
write.tree(bifurcating_tree_lb_10_change,file = "/project/exaptation/fast_approx_pruning/test_data/test_lb_10_change/test_lb_10_change_bi.fastree")
