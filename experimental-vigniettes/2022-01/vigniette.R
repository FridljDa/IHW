#http://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
library("DESeq2")
library("dplyr")
data("airway", package = "airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex) %>% DESeq
deRes <- as.data.frame(results(dds))

colnames(deRes)

library("IHW")
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes, alpha = 0.1)

rejections(ihwRes)

head(adj_pvalues(ihwRes))
sum(adj_pvalues(ihwRes) <= 0.1, na.rm = TRUE) == rejections(ihwRes)

padjBH <- p.adjust(deRes$pvalue, method = "BH")
sum(padjBH <= 0.1, na.rm = TRUE)

head(weights(ihwRes))

c(nbins(ihwRes), nfolds(ihwRes))

weights(ihwRes, levels_only = TRUE)

plot(ihwRes)
plot(ihwRes, what = "decisionboundary") 

library("ggplot2")
gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg

gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)

ihwResDf <- as.data.frame(ihwRes)
colnames(ihwResDf)

ihwBonferroni <- ihw(pvalue ~ baseMean, data = deRes, alpha = 0.1, adjustment_type = "bonferroni")
weights(ihwBonferroni)

#3 diagnostics
deRes <- na.omit(deRes)
deRes$geneid <- as.numeric(gsub("ENSG[+]*", "", rownames(deRes)))

# set up data frame for ggplotting
rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes), 
                 covariate_type="base mean"),
      data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$geneid)/nrow(deRes), 
                 covariate_type="gene id")) %>%
  ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + 
  facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)

deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 8)

ggplot(deRes, aes(x=pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup, nrow = 2)

ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step") 

deRes$lfcGroup <- groups_by_filter(abs(deRes$log2FoldChange),8)
ggplot(deRes, aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ lfcGroup, nrow=2)

ggplot(deRes, aes(x = pvalue, col = lfcGroup)) + stat_ecdf(geom = "step") 
