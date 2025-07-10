xlabAngle = 0
xlabhjust = 0.5
xlabvjust = 0.5
ylabAngle = 0
ylabhjust = 0.5
ylabvjust = 0.5
axisLabSize = 16
titleLabSize = 16
subtitleLabSize = 12
captionLabSize = 12
legendPosition = 'top'
legendLabSize = 10
legendIconSize = 3.0

theme_set(
  theme_bw(base_size=24) +
    theme(
      legend.background=element_rect(),
      
      plot.title=element_text(angle=0, size=16,
                              face='bold', vjust=1),
      plot.subtitle=element_text(angle = 0, size = 12,
                                 face = 'plain', vjust = 1),
      plot.caption=element_text(angle = 0, size = 12,
                                face = 'plain', vjust = 1),
      
      axis.line = element_line(size=1.5, colour = 'black'),
      
      axis.text.x=element_text(angle = xlabAngle, size = 16,
                               hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y=element_text(angle = ylabAngle, size = 16,
                               hjust = ylabhjust, vjust = ylabvjust),
      axis.title=element_text(size=16),
      
      legend.position=legendPosition,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.key=element_blank(),
      legend.key.size=unit(0.5, 'cm'),
      legend.text=element_text(size=legendLabSize),
      
      title=element_text(size=legendLabSize),
      legend.title=element_blank())
)