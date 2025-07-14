xlabAngle = 0
xlabhjust = 0.5
xlabvjust = 0.5
ylabAngle = 0
ylabhjust = 0.5
ylabvjust = 0.5
axisLabSize = 10
titleLabSize = 12
subtitleLabSize = 10
captionLabSize = 10
legendPosition = 'top'
legendLabSize = 10
legendIconSize = 3.0

theme_set(
  theme_bw(base_size=24) +
    theme(
      legend.background=element_rect(),
      
      plot.title=element_text(angle=0, size=titleLabSize,
                              face='bold', vjust=1),
      plot.subtitle=element_text(angle = 0, size = subtitleLabSize,
                                 face = 'plain', vjust = 1),
      plot.caption=element_text(angle = 0, size = captionLabSize,
                                face = 'plain', vjust = 1),
      
      axis.line = element_line(size=1.5, colour = 'black'),
      
      axis.text.x=element_text(angle = xlabAngle, size = axisLabSize,
                               hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y=element_text(angle = ylabAngle, size = axisLabSize,
                               hjust = ylabhjust, vjust = ylabvjust),
      axis.title=element_text(size=titleLabSize),
      
      legend.position=legendPosition,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.key=element_blank(),
      legend.key.size=unit(0.5, 'cm'),
      legend.text=element_text(size=legendLabSize),
      
      title=element_text(size=legendLabSize),
      legend.title=element_blank())
)