library(ggpubr)
library(ggplot2)

load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/script5a.RData")


Kmeans_silhouette_all_edition_ada2<- Kmeans_silhouette
Kmeans_silhouette_all_edition_ada2$k<-as.character(Kmeans_silhouette_all_edition_ada2$k)
# lollipop plot

title_size=40
axis_size=15

a1<-ggplot(Kmeans_silhouette_all_edition_ada2, aes(x=k, y=mean_silhouette_score)) +
  geom_segment( aes(x=fct_inorder(k), xend=k, y=0, yend=mean_silhouette_score)) +geom_point()  +
  geom_point(data = Kmeans_silhouette_all_edition_ada2[which.max(Kmeans_silhouette_all_edition_ada2$mean_silhouette_score), ], color="red")+
  ylab("Mean Silhouette Score")+ggtitle("a")+theme(plot.title = element_text(size = title_size, face = "bold"),axis.text=element_text(size=axis_size),
                                                   axis.title=element_text(size=axis_size,face="bold"))

rm(list=setdiff(ls(), c("Kmeans_silhouette_all_edition_ada2","a1")))


load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/script5a_5thed.RData")

Kmeans_silhouette_5th_edition_ada2<- Kmeans_silhouette
Kmeans_silhouette_5th_edition_ada2$k<-as.character(Kmeans_silhouette_5th_edition_ada2$k)

rm(list=setdiff(ls(), c("Kmeans_silhouette_all_edition_ada2","a1","Kmeans_silhouette_5th_edition_ada2")))

a2<-ggplot(Kmeans_silhouette_5th_edition_ada2, aes(x=k, y=mean_silhouette_score)) +
  geom_segment( aes(x=fct_inorder(k), xend=k, y=0, yend=mean_silhouette_score)) +geom_point()  +
  geom_point(data = Kmeans_silhouette_5th_edition_ada2[which.max(Kmeans_silhouette_5th_edition_ada2$mean_silhouette_score), ], color="red")+
  ylab("Mean Silhouette Score")+ggtitle("b")+theme(plot.title = element_text(size = title_size, face = "bold"),axis.text=element_text(size=axis_size),
                                                    axis.title=element_text(size=axis_size,face="bold"))

load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/script-5B.RData")

Kmeans_silhouette_all_edition_v3<- Kmeans_silhouette
Kmeans_silhouette_all_edition_v3$k<-as.character(Kmeans_silhouette_all_edition_v3$k)
rm(list=setdiff(ls(), c("Kmeans_silhouette_all_edition_ada2","a1","Kmeans_silhouette_5th_edition_ada2","a2",
                        "Kmeans_silhouette_all_edition_v3")))

a3<-ggplot(Kmeans_silhouette_all_edition_v3, aes(x=k, y=mean_silhouette_score)) +
  geom_segment( aes(x=fct_inorder(k), xend=k, y=0, yend=mean_silhouette_score)) +geom_point()  +
  geom_point(data = Kmeans_silhouette_all_edition_v3[which.max(Kmeans_silhouette_all_edition_v3$mean_silhouette_score), ], color="red")+
  ylab("Mean Silhouette Score")+ggtitle("c")+theme(plot.title = element_text(size = title_size, face = "bold"),axis.text=element_text(size=axis_size),
                                                    axis.title=element_text(size=axis_size,face="bold"))

load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/script-5b5thed-aug21.RData")
Kmeans_silhouette_5th_edition_v3<- Kmeans_silhouette
Kmeans_silhouette_5th_edition_v3$k<-as.character(Kmeans_silhouette_5th_edition_v3$k)


rm(list=setdiff(ls(), c("Kmeans_silhouette_all_edition_ada2","a1","Kmeans_silhouette_5th_edition_ada2","a2","a3","Kmeans_silhouette_5th_edition_v3",
                        "Kmeans_silhouette_all_edition_v3")))
a4<-ggplot(Kmeans_silhouette_5th_edition_v3, aes(x=k, y=mean_silhouette_score)) +
  geom_segment( aes(x=fct_inorder(k), xend=k, y=0, yend=mean_silhouette_score)) +geom_point()  +
  geom_point(data = Kmeans_silhouette_5th_edition_v3[which.max(Kmeans_silhouette_5th_edition_v3$mean_silhouette_score), ], color="red")+
  ylab("Mean Silhouette Score")+ggtitle("d")+theme(plot.title = element_text(size = title_size, face = "bold"),axis.text=element_text(size=axis_size),
                                                   axis.title=element_text(size=axis_size,face="bold"))


p5<-ggarrange(a1, a2,a3,a4,nrow = 2,ncol = 2)
ggsave(p5, filename ="Paper/plots/Kmeans_Silhouette_Score_Embeddings_all.png",, height = 50, width = 120, units = "cm")
