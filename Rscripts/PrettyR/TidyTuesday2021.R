#=================================================
#2021-01-05 unfinished 
#=================================================
#TODO
library(tidyverse)
library(textreadr)
library(ggforce)
library(glue)
library(systemfonts)
library(ggtext)
library(pdftools)

theme_set(theme_light(base_family = "Roboto Condensed", base_size = 20))
theme_update(
  axis.title = element_blank(),
  axis.text.x = element_text(color = "grey70", size = 18, face = "bold", hjust = .15),
  axis.text.y = element_markdown(margin = margin(r = 10)), 
  axis.ticks.x.top = element_line(
    size = .6,
    arrow = arrow(length = unit(.4, "lines"), ends = "first", type = "closed")
  ),
  axis.ticks.x.bottom = element_line(
    size = .6,
    arrow = arrow(length = unit(.4, "lines"), ends = "last", type = "closed")
  ),
  axis.ticks.y = element_blank(),
  axis.ticks.length.x = unit(1.3, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(fill = NA, color = NA),
  legend.position = "none",
  plot.caption = element_markdown(color = "grey70", size = 14, lineheight = 1.15,
                                  margin = margin(t = 20)),
  plot.margin = margin(30, 45, 12, 45)
)



#=================================================
#2019-11-12 unfinished
#=================================================
#TODO
library(tidyverse)
library(igraph)
library(ggraph)

setwd('D:\\R\\Rproject\\scripts\\tidytuesday')
data <- read.csv('loc_cran_packages.csv')

# most popular programming languages from TIOBE Index (Nov. 2019) found in data
# (only languages with position <= 16 are considered)
popular_languages <- c(
  'Java', 'C', 'Python', 'C++', 'C#', 'Visual Basic', 'JavaScript', 'PHP', 'SQL', 'Ruby', 'Objective C++', 'Assembly', 'R'
)


# number of packages to display
number_of_pkgs <- 300

# find largest packages written in popular languages
top_packages <- data %>% filter(language %in% popular_languages) %>% 
  group_by(pkg_name) %>% summarise(total_code = sum(code)) %>% 
  arrange(desc(total_code)) %>% head(number_of_pkgs) %>% select(pkg_name, total_code)

# all popular languages per package
data %>%filter(pkg_name %in% top_packages$pkg_name,language %in% popular_languages) %>% 
  arrange(pkg_name, desc(code))

# all popular languages per package
top_languages_per_pkg <- data %>%
  filter(
    pkg_name %in% top_packages$pkg_name,
    language %in% popular_languages
  ) %>%
  arrange(pkg_name, desc(code)) %>%
  group_by(pkg_name) %>%
  mutate(
    main = row_number() == 1, # main language of package should be opaque
    total_code = sum(code)
  ) %>%
  ungroup() %>%
  select(language, pkg_name, code, total_code, main)

# only following languages found in given packages
# Add parentheses to assign the result and print the result at the same time
(top_languages <- top_languages_per_pkg %>%
    pull(language) %>%
    unique %>%
    sort)

top_language_colors <- c(
  '#efb306',
  '#eb990c',
  '#e8351e',
  '#cd023d',
  '#852f88',
  '#4e54ac',
  '#0f8096',
  '#7db954',
  '#17a769',
  '#000000'
)

names(top_language_colors) <- c(
  'Assembly',
  'C',
  'C++',
  'JavaScript',
  'Java',
  'R',
  'Python',
  'Ruby',
  'SQL',
  'All'
)

edges1 <- top_languages_per_pkg %>%
  transmute(from = language, to = pkg_name, total_code = code, main)

edges2 <- top_languages_per_pkg %>%
  count(language, wt = code, name = 'total_code') %>%
  transmute(
    from = '',
    to = language,
    total_code,
    main = TRUE
  )

vertices1 <- top_languages_per_pkg %>%
  filter(main) %>%
  transmute(
    node = pkg_name, language, total_code, level = 1
  )
vertices2 <- edges2 %>%
  transmute(
    node = to, language = to, total_code, level = 2
  )

vertices3 <- tibble(
  node = '', language = NA, total_code = 0, level = 3
)

vertices <- bind_rows(vertices1, vertices2, vertices3) %>%
  mutate(
    radius = total_code**(1.8), # scaling circles
    language = factor(language, names(top_language_colors))
  ) %>%
  arrange(level, language, node)

graph <- graph_from_data_frame(edges, vertices = vertices)

# create custom layout by updating existing circle layout
layout <- create_layout(graph, layout = 'circle')


#=================================================
#[2021-03-23]
#@Finished time:2021-03-31
#@Data description:https://github.com/rfordatascience/tidytuesday/tree/master/data/2021/2021-03-23
#=================================================
library(tidyverse)
library(tidytuesdayR)
library(tidytext)
library(ggtext)
library(lubridate)
library(ggthemes)

un_votes <- read.csv("unvotes.csv",sep=",")
issues <- read.csv("issues.csv",sep=",")
roll_calls <- read.csv("roll_calls.csv",sep=",")

#------------------------柱状图+线图------------------------
roll_calls <- roll_calls %>% 
  mutate(year = year(date)) %>%    #year(date) 将1946-01-01转换为1946
  left_join(issues)

roll_calls %>%
  group_by(year) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = year, y = n)) +
  geom_line()

un <- roll_calls %>%
  filter(!is.na(issue)) %>%   #将issue一列为空的去掉
  count(year, issue, sort = T) %>% arrange(year)
un_plot <- ggplot(un,aes(x=year,y=n)) + geom_col(aes(fill=issue)) +
  geom_line(aes(x=year,y=n),color='black',data=roll_calls %>% count(year)) +
  scale_fill_ptol() + #换种颜色，来自ggthems包
  theme_bw()+theme(legend.title = element_blank(),legend.position='top')

ggsave("un_plot.pdf", un_plot, height = 8, width = 12, units = "in", dpi = 300)

#-----------------------柱状图-----------------------------
colors <-c("#E41A1C","#1E90FF","#FF8C00","#4DAF4A","#984EA3",
           "#40E0D0","#FFC0CB","#00BFFF","#FFDEAD","#90EE90",
           "#EE82EE","#00FFFF","#F0A3FF", "#0075DC", 
           "#993F00","#4C005C","#2BCE48","#FFCC99",
           "#808080","#94FFB5","#8F7C00","#9DCC00",
           "#C20088","#003380","#FFA405","#FFA8BB",
           "#426600","#FF0010","#5EF1F2","#00998F",
           "#740AFF","#990000","#FFFF00")

bar <- roll_calls %>%
  filter(!is.na(issue)) %>%
  count(year, issue, sort = T) %>%
  arrange(year) %>%
  ggplot(aes(x = year, y = n)) +
  geom_col(aes(fill = issue))+
  scale_fill_manual(values = colors)+labs(fill="")+
  labs(x=NULL,y=NULL)+
  theme_bw()+theme(legend.position = "top")

ggsave("bar.pdf", bar, height = 8, width = 12, units = "in", dpi = 300)



#=================================================
#[2021-03-02]
#@Finished time:2021-04-03
#@Data description:https://github.com/rfordatascience/tidytuesday/tree/master/data/2021/2021-03-02
#=================================================
library(tidyverse)
library(scales)
library(patchwork)
library(ggsci)

tues_data <- tidytuesdayR::tt_load(2021, week = 10)

youtube <- tues_data$youtube

df_longer <- youtube %>% 
  pivot_longer(cols = funny:use_sex, names_to = "trait",    #pivot_longer 长数据转变为宽数据[tidyr]
               values_to = "present", 
               values_ptypes = list("present" = integer())) %>% 
  select(year, brand, trait, present)

heatmap <- df_longer %>% 
  group_by(year, trait) %>% 
  summarise(total = sum(present)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mutate(prop = total / sum(total)) %>% 
  ggplot(aes(year)) +
  geom_tile(aes(y = trait, fill = prop),   #geom_tile 绘制热图
            color = "white", size = 0.5) + #size 指定热图之间的分割线的粗细
  scale_x_continuous(expand = c(0, 0)) +   # 删去图片在x轴左右的留白
  scale_fill_material("cyan", guide = guide_colorsteps(position = "bottom"), 
                      labels = percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(
    plot.margin = margin(0, 20, 20, 20),  #调整图片在页面显示的位置
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(2, "cm"), #设置图例的长短粗细
    axis.title = element_blank(),
    axis.text.x = element_text(size = 15),  #改变x轴字体大小
    axis.text.y = element_text(size = 15),    #改变y轴字体大小
    plot.caption = element_text(size = 13, hjust = 0.5)
  )

bar_plot <-
  df_longer %>% 
  group_by(year) %>% 
  summarise(total = sum(present)) %>% 
  ggplot(aes(year, total)) +
  geom_col(fill = "#B2EBF2FF") +
  geom_rect(aes(xmin = year - 0.45, xmax = year + 0.45, 
                ymin = total - 5, ymax = total), fill = "#00838EFF") +  #geom_rect 添加阴影区域
  geom_text(aes(y = total - 2.5, label = total), 
            color = "#DFF7F9FF",size = 6.5, fontface = "bold") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Evolution of Superbowl Commercials",
    subtitle = "The bars represent the number of ads per year between 2000-2020.\nThe heatmap shows the distribution of themes prevalent in the commercials."
  ) +
  theme_void() +  #theme_void去除x轴y轴的字体显示
  theme(
    plot.title = element_text(size = 30, hjust = 0.5, margin = margin(b = 10)),  #修改主标题
    plot.subtitle = element_text(size =15, hjust = 0.5, margin = margin(b = 15)), #修改次级标题
    plot.margin = margin(t = 20)
  )

result_plot <- bar_plot/heatmap

ggsave("20210302_result.pdf", result_plot, height = 12, width = 16, units = "in", dpi = 300)
