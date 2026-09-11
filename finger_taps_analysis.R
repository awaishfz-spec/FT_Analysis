library(signal)
library(tidyverse)
library(zoo)
library(pracma)
library(optotools)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


## Functions -------
# Optotrak Files
read_and_filter <- function(file){
  
  movement <- read_ndi(file) 
  dt <- 1/ movement$freq
  
  df <- as.data.frame(movement) %>% 
    mutate_all(na.approx, na.rm = FALSE) %>% # lerp
    mutate(frame = 0:(movement$numframes-1))
  
  # pad data
  paddata <- 100
  
  data <- tibble(frame = (-paddata):(movement$numframes+paddata-1))
  
  data <- data %>%
    left_join(y = df, by = "frame") %>%
    fill(-frame, .direction = "down") %>% # padded NAs replaced with first/last position
    fill(-frame, .direction = "up")
  
  nyq_f <- movement$freq/2
  
  bf_low <- butter(4, 5/nyq_f , type="low",plane = 'z')
  
  # filter  
  data <- data %>% 
    mutate_at(vars(-frame), function(x) filtfilt(filt = bf_low, x = x)) %>% 
    dplyr::filter(frame >= 0, frame < movement$numframes) %>% 
    mutate(time = frame * dt) %>% 
    select(-frame)
  
  return(data)
}

read_and_filter_big <- function(file){
  
  movement <- read_csv(file) 
  dt <- 1/ 60
  numframes = length(movement$m0_x)
  df <- as.data.frame(movement) %>% 
    mutate_at(vars(m0_x:m1_z), funs(ifelse(. == -999, NA, .))) %>% 
    mutate_at(vars(m0_x:m1_z),na.approx, na.rm = FALSE) %>% # lerp
    mutate(frame = 0:(numframes-1))
  # pad data
  paddata <- 100
  
  data <- tibble(frame = (-paddata):(numframes+paddata-1))
  
  data <- data %>%
    left_join(y = df, by = "frame") %>%
    fill(-frame, .direction = "down") %>% # padded NAs replaced with first/last position
    fill(-frame, .direction = "up")
  
  nyq_f <- 30
  
  bf_low <- butter(4, 5/nyq_f , type="low",plane = 'z')
  
  # filter  
  data <- data %>% 
    mutate_at(vars(-frame), function(x) filtfilt(filt = bf_low, x = x)) %>% 
    dplyr::filter(frame >= 0, frame < numframes) %>% 
    mutate(time = frame * dt) %>% 
    select(-frame)
  
  return(data)
}


#Bigkat Files
read_and_time <- function(file){
  
  movement <- read_csv(file) 
  dt <- 1/ 60
  numframes = length(movement$m0_x)
  df <- as.data.frame(movement) %>% 
    mutate(frame = 0:(numframes-1))
  
  data <- df %>% 
    mutate(time = frame * dt) %>% 
    select(-frame)
  
  return(data)
}


# Sampling Time
delta_t <- 1/60

# raincloud 
raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))

## Import data and extract relevant parameters ---------

big_file_filt <- tibble(files = list.files("BIGKAT_csv", full.names= TRUE, recursive = TRUE, pattern = "^3d_unfil(.*)csv$"))

bigkat_filt <- big_file_filt %>% 
  mutate(data = map(files, read_and_filter_big)) %>% 
  unnest(data, .drop = FALSE) %>% 
  separate(files, c("device", "ppid", "records","file_num"), sep = "/",extra='drop') %>% 
  #select(-records) %>% 
  separate(device, c("device","csv"),sep = "_") %>%
  select(-csv,-ppid) %>%
  separate(records, c("ppid","hand"),sep = "_") %>% 
  group_by(ppid, hand) %>% 
  mutate(trial_num = as.integer(as.factor(file_num))) %>% 
  select(-file_num) %>% 
  rename_(.dots=setNames(names(.), tolower(gsub("m0", "thumb", names(.))))) %>% 
  rename_(.dots=setNames(names(.), tolower(gsub("m1", "finger", names(.))))) %>% 
  select(device,ppid,hand,trial_num,everything())%>% 
  data.frame()

#bigkat_filt %>% group_by(ppid, hand, trial_num) %>% summarise() 

library(data.table)
bigkat_filt <- setnames(bigkat_filt, old = c('thumb_x','thumb_y','finger_x','finger_y'), new = c('thumb_y','thumb_x','finger_y','finger_x'))
bigkat_filt <- bigkat_filt %>%
  mutate(thumb_y = -1 * thumb_y,
         finger_y = -1 * finger_y,
         thumb_x = -1 * thumb_x,
         finger_x = -1 * finger_x)



opto_files <- tibble(files = list.files("OPTOTRAK",
                                        full.names = TRUE,
                                        recursive = TRUE, pattern = "^C#(.*)P70|P71|P72|P74|C01|C02|C03|P75|P77|P78|P82|P83|P89|
                                        |C05|C06|C07|C08|C09|C10|C11|C12|C13"))

optotrak <- opto_files %>% 
  mutate(move_data = map(files, read_and_filter), device = "OPTOTRAK") %>% 
  mutate(file_num = substr(files, 12, 18)) %>% 
  left_join(read_csv("file_number_updated_cor.csv"), by = "file_num") %>% 
  unnest(move_data, .drop = FALSE) %>%
  rename_(.dots=setNames(names(.), tolower(gsub("m1", "thumb", names(.))))) %>% 
  rename_(.dots=setNames(names(.), tolower(gsub("m2", "finger", names(.))))) %>% 
  group_by(ppid, hand,device) %>% 
  mutate(trial_num = as.integer(as.factor(file_num))) %>% 
  #mutate(trial_num =  file_number - min(file_number)) %>% 
  #mutate(trial_num=as.integer(trial_num)) %>% 
  select(-files,-file_num,-file_number,-ext) %>% 
  select(device,ppid,hand,trial_num,everything())

optotrak<-optotrak[complete.cases(optotrak),]

#optotrak <- setnames(optotrak, old = c('thumb_x','thumb_y','finger_x','finger_y'), new = c('finger_y','finger_x','thumb_y','thumb_x'))

# had to delete after 17 sec for P65, for left hand because data was collected for 70 secs by optotrak
optotrak_P65_left <- optotrak %>% 
  dplyr::filter(ppid == 'P65') %>% 
  dplyr::filter(hand == 'Left') %>% 
  dplyr::filter(time < 17)
#Get data for P65 & for right hand
optotrak_P65_right <- optotrak%>% 
  dplyr::filter(ppid == 'P65') %>% 
  dplyr::filter(hand == 'Right')
# combine data for P65  for both hands
Optotrak_P65 <- bind_rows(optotrak_P65_right,optotrak_P65_left)

# Optotrak data without P65
optotrak_excep_P65 <- optotrak %>% 
  dplyr::filter(ppid != 'P65')
# combine whole optotrak data
optotrak <- bind_rows(optotrak_excep_P65,Optotrak_P65)

Fingertapping_raw <- bind_rows(bigkat_filt, optotrak) 
#write_csv('Fingertapping_raw.csv')


## For matlab
df <- Fingertapping_raw
df$fileName <- paste(df$device, "_", df$ppid, "_", df$hand, "_", df$trial_num) 

df <- df %>%
  select(-device,-ppid,-hand,-trial_num,-time) %>% 
  select(fileName,everything())

##save file in csvfiles_big_op
lapply(split(df, df$fileName),
       function(x) write.csv(x, file = paste0(x$fileName[1], ".csv",sep=""),row.names = FALSE))



## Processing-------
Fingertapping_raw <- read_csv("Fingertapping_raw.csv")

#delete first second of the record
FT_raw <- Fingertapping_raw  %>% dplyr::filter(time >= 1) %>% 
  group_by(ppid,hand,trial_num,device) %>% 
  mutate(t = (row_number() * delta_t)-delta_t,
         distance_thumb_finger = ((finger_x - thumb_x)^2 + (finger_y - thumb_y)^2 + (finger_z - thumb_z)^2)^0.5,
         speed_finger = abs(c(NA, sqrt((diff(finger_x)^2 + diff(finger_y)^2 + diff(finger_z)^2))/diff(t)))) %>% 
        #finger_peak_speed = max(speed_finger, na.rm = TRUE),
  select(-time,-starts_with('th'),-starts_with('fi'))
  
fingertapping_movement <- FT_raw  %>%
  group_by(ppid,hand,trial_num,device) %>% 
  mutate(amp_sig = distance_thumb_finger / max(distance_thumb_finger),
         speed_thumb_finger = abs(c(NA ,diff(amp_sig)/diff(t))),
         ave_speed = mean(speed_thumb_finger,na.rm = TRUE)) 

ft_summary <- fingertapping_movement  %>%
  group_by(ppid,hand,trial_num,device) %>% 
  summarize(amp_sig = mean(amp_sig),
         ,
         ave_speed = mean(ave_speed,na.rm = TRUE),
         distance_thumb_finger= mean(distance_thumb_finger),
         speed_finger = mean(speed_finger,na.rm = TRUE))


amplitude_CV <- function(df_subset){
  amp_signal <- df_subset$amp_sig
  amp_signal <- amp_signal - mean(amp_signal)
  win_len <- 60
  energy_list <- 0

  l <- length(amp_signal)
  ra <- l - win_len

  for (i in 0:ra){

    temp <- amp_signal[i:win_len+i]
    energy_list[i] <- max(temp,na.rm = TRUE) - min(temp,na.rm = TRUE)
  }

  df <- tibble(energy_CV = sd(energy_list)/mean(energy_list))

  return(df)
}

ftm <- fingertapping_movement  %>%
  group_by(ppid,hand,trial_num,device) %>%
  do(amplitude_CV(.))

ftm <- ftm %>%  mutate( part = ifelse(grepl('C',ppid), 'Control', "Patient"))
ftm_big <- ftm %>% filter(device == 'BIGKAT')


ggplot(ftm, aes(x = part, y = energy_CV, fill = device))+
  #geom_line()+
  geom_boxplot(outlier.shape  = NA)+
  #geom_point(size=0.7, position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1),show.legend = FALSE)+
  facet_grid(~hand)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black'))+
  theme(axis.ticks = element_blank())+
  scale_fill_manual(values = c("tomato", "green"))+
  labs(x = "", y = "Amplitude (mm)", fill = "Device",tag = 'a')+
  theme(legend.text=element_text(size=7))+
  theme(legend.title=element_text(size=9)) 




## FFT of amplitude

ampp <- fingertapping_movement %>% ungroup() %>% filter(device=='BIGKAT') %>% 
  filter(ppid=='C01') %>%
  filter(hand=='Left') %>%
  filter(trial_num=='1') %>% 
  select(amp_sig)



  
ftm <- fingertapping_movement


for (i in 0:length(ftm$amp_sig - 60)){
  temp = 
}

amp_win <- ftm %>% 
  group_by(ppid,hand,trial_num,device) %>%
  lapply()

amp_win <- ftm %>% 
  group_by(ppid,hand,trial_num,device) %>%
  sapply(split(ftm$amp_sig, ceiling(seq_along(ftm$amp_sig)/16)) %>% mean)

max_amp <- split(fingertapping_movement$amp_sig, ceiling(seq_along(fingertapping_movement$amp_sig)/16))
              
max_amp <-sapply(split(fingertapping_movement$amp_sig, ceiling(seq_along(fingertapping_movement$amp_sig)/16)),mean )
window_length = 60;
for (i = 0 : length(fingertapping_movement$signal) - window_length):
  temp = signal(i : window_length+i-1),
  
  energy_list(i) = max(temp) - min(temp)
end


# graphs

ggplot(fingertapping_movement %>% dplyr::filter(ppid =='P77')%>% dplyr::filter(hand =='Right')%>% dplyr::filter(trial_num =='1'),aes(x=time,y=abs(distance_thumb_finger),color = device))+
  geom_line()+
  #scale_x_continuous(limits = c(6,12))+
  #facet_wrap(hand~trial_num)+
  theme_bw()

ggplot(Fingertapping_raw %>% dplyr::filter(ppid =='P64') %>% dplyr::filter(device =='BIGKAT') , aes(x = time, y = finger_x), color = 'black') +
  geom_line() + 
  geom_line(aes(y=thumb_x),color='red')+
  facet_wrap(hand~ trial_num, scales = "free") + 
  #scale_x_continuous(limits = c(6,10))+ 
  theme_bw()+
  labs(x="Time (s)",y = "Y-axis")  

ggplot(Fingertapping_raw %>% dplyr::filter(ppid =='P82') , aes(x = time, y = finger_x), color = 'green') +
  geom_line() + 
  geom_line(aes(y=thumb_x),color='red')+
  facet_wrap(hand*device~ trial_num, scales = "free") + 
  #scale_x_continuous(limits = c(6,10))+ 
  theme_bw()+
  labs(x="Time (s)",y = "Y-axis")  

ggplot(fingertapping_movement%>% dplyr::filter(ppid =='P64') , aes(x = ppid, y = abs(distance_thumb_finger), color = device)) +
  #geom_boxplot()+ 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .5, position = position_dodge(width = 1))+
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(width = 1))+
  #geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1),show.legend = FALSE, alpha = .3)+
  theme_bw()+
  facet_wrap(trial_num~ppid)+
  labs(x='Participants',y = "Max Amplitude (mm)") 
#coord_cartesian(ylim = c(0, 200))
#+ 
#  stat_summary(aes(group = ppid), fun.y = mean, geom = "point", size = 2))
ggsave("max_amplitude.png", graph, type = "cairo-png")

ggplot(fingertapping_movement %>% group_by(hand) %>% dplyr::filter(ppid =='P64'), aes(x = ppid, y = velocity_thumb_finger, color = device)) +
  #geom_boxplot()+ 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .5, position = position_dodge(width = 1))+
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(width = 1))+
  #geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1),show.legend = FALSE, alpha = .3)+
  theme_bw()+
  facet_wrap(trial_num~ppid)



# Processing



df_summary <- fingertapping_movement %>%
  group_by(ppid, trial_num,device) %>% 
  summarise(area_under_curve = sum(delta_t * rollmean(distance_thumb_finger,2 ))) %>% 
  ungroup()


df_fft <- fingertapping_movement %>%
  select(ppid, trial_num, time,distance_thumb_finger) %>%   
  dplyr::filter(time > 1.3, time < 12) %>% 
  group_by(ppid, trial_num) %>% 
  mutate(pow= Mod(fft(distance_thumb_finger)),freq = delta_t * row_number()/n()) %>% 
  ungroup() 

ggplot(df_fft, aes(pow))+
  geom_histogram(bins = 10)+
  facet_grid(ppid~.)

ggplot(df_fft , aes(x=freq, y= pow, color=ppid,fill=ppid))+
  geom_line()+
  facet_grid(~ppid,scales='free')+
  stat_smooth(geom = 'area',
              method = 'loess',
              span = .9,
              alpha = .2)+ #+
  #scale_x_continuous(limits = c(4, 22))+
  scale_y_continuous(limits = c(0, 1800))

df_status <- df_fft %>% 
  mutate(status = 'TRUE')


df_long <- fingertapping_movement %>%
  ungroup() %>% 
  select(device,ppid, hand, trial_num, time, distance_thumb_finger,starts_with("fin"),starts_with("th")) %>% 
  select(-finger_speed_vec,-finger_peak_speed,-finger_peak_speed_time,-thumb_speed_vec,-thumb_peak_speed,-thumb_peak_speed_time) %>% 
  gather(variable, value,starts_with("fin"),starts_with("th")) %>% 
  #mutate(variable = gsub("_mov", "", variable)) %>% 
  separate(variable, c("part", "axis"), sep = "_") %>% 
  spread(axis, value) 


ggplot(df_long %>% group_by(ppid) , aes(x = time, y =- y, color = device)) +
  geom_line() + 
  #facet_wrap(~ trial_num, scales = "free") + 
  scale_x_continuous(limits = c(6, 9))+ 
  theme_bw()+
  facet_wrap(trial_num ~ hand)
labs(x="Time (s)",y = "Y-axis")  

df_fft_part <- df_long %>% 
  group_by(ppid, trial_num,part) %>% 
  dplyr::filter(time > 1.3, time < 12) %>% 
  mutate(pow_y= Mod(fft(y)) ,freq = delta_t * row_number()/n()) %>% 
  mutate(pow_th_fi= Mod(fft(distance_thumb_finger)) ,freq = delta_t * row_number()/n()) %>% 
  ungroup() 

ggplot(df_fft_part, aes(pow_y))+
  geom_histogram()+
  facet_grid(ppid~.)

ggplot(df_fft_part, aes(x=freq, y= pow_y, color=part,fill=part))+
  geom_line()+
  facet_grid(part~ppid,scales='free')+
  stat_smooth(geom = 'area',
              method = 'loess',
              span = 11.5,
              alpha = .2)+ #+
  scale_x_continuous(limits = c(0, 60))+
  scale_y_continuous(limits = c(0, 1000))

ggplot(df_fft_part, aes(x=freq, y= pow_th_fi, color=ppid,fill=ppid))+
  geom_line()+
  facet_grid(~ppid,scales='free')+
  stat_smooth(geom = 'area',
              method = 'loess',
              span = 11.5,
              alpha = .2)+ #+
  scale_x_continuous(limits = c(0, 30))+
  scale_y_continuous(limits = c(0, 1000))


set_peaks_column <- function(x){
  peaks <- findpeaks(x, nups = 5, ndowns = 5)
  peaks_idx <- peaks[,2]
  
  is_peak <- rep(FALSE, length(x))
  
  is_peak[peaks_idx] = TRUE
  return(is_peak)
}


df_peaks <- fingertapping_movement %>% 
  select(device,ppid, hand,trial_num, time, distance_thumb_finger) %>% 
  group_by(device,hand,ppid, trial_num) %>% 
  mutate(velocity_thumb_finger = c(NA, diff(distance_thumb_finger) / delta_t)) %>% 
  dplyr::filter(!is.na(velocity_thumb_finger)) %>% 
  mutate(is_peak = set_peaks_column(velocity_thumb_finger),
         is_trough = set_peaks_column(-velocity_thumb_finger)) %>% 
  dplyr::filter(is_peak | is_trough) %>% 
  mutate(feature = ifelse(is_peak, "peak", "trough"),
         half_time_period = c(NA, diff(time)),
         Time_Period = 2* half_time_period,
         freq = 1 / Time_Period) #%>% 
select(-is_peak, -is_trough)


df_freq <- df_peaks %>% 
  select(ppid,trial_num,time, feature) %>% 
  group_by(ppid,trial_num) %>% 
  mutate(Time_Period = c(NA,NA, diff(time,lag=2)),
         freq = 1 / Time_Period) %>% 
  ungroup()

ggplot(df_peaks %>% group_by(feature),aes(x=feature, y=freq,color=feature))+
  geom_boxplot(position = position_dodge(width = 0.5))+
  facet_wrap(device~ppid*hand)


ggplot(df_freq ,aes(freq,color=ppid,fill=ppid))+
  geom_histogram()+
  facet_grid(ppid~.)

ggplot(df_freq %>% group_by(ppid,trial_num) ,aes(x=time, y=freq,color=ppid))+
  geom_line()+
  #stat_summary(fun.data = mean_se, geom = "point", size = 2)+
  facet_grid(ppid~.)+
  theme_bw()


ggplot(df_freq %>% group_by(ppid,trial_num) ,aes(x=time, y=freq,color=ppid))+
  geom_boxplot()+
  stat_summary(fun.data = mean_se, geom = "point", size = 2)+
  facet_grid(~ppid)+
  theme_bw()


ggplot(df_freq ,aes(x=time, y=freq,color=ppid,group(ppid,trial_num)))+
  geom_boxplot()+
  stat_summary(aes(group = trial_num), fun.data = mean_se, geom = "point", size = 2)+
  facet_grid(~ppid)+
  theme_bw()
stat_summary(aes(group = trial_num), fun.data = mean_se, geom = "point", size = 2)

#error bar

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

sumld1 <- df_freq %>% 
  #select(-time) %>% 
  dplyr::filter(!is.na(freq)) %>% 
  group_by(ppid,trial_num,feature) %>% 
  summarise_at(vars(freq), funs(median = median, lower = lb, upper = ub, freq = mean))

sumld1

ggplot(df_freq %>% dplyr::filter(!is.na(freq)), aes(x = feature, y = freq, fill = feature,group=feature)) + 
  geom_flat_violin(position = position_nudge(x = .1), alpha = .8) + 
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  #stat_summary(fun.data = mean_se, geom = "point", size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper,group=feature),
                data = sumld1,
                position = position_nudge(x = -0.2),
                width = 0.1) +
  #geom_point(data = sumld1, position = position_nudge(x = -0.4), size = .8, alpha = .8) +
  #geom_point(position = position_jitter(width = .1), size = .25) +
  #scale_x_continuous(breaks = c(), name = "") +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  facet_grid(~ppid, switch = 'x') + 
  theme_bw() +
  theme(strip.background = element_rect(color = 'black',fill = 'NA')) +
  raincloud_theme+
  #theme(panel.border = element_rect(fill=NA))+
  xlab('')+
  ylab('Error (mm)')+
  ggsave('Graph21.png', type = "cairo", height = 5, width = 5, dpi = 600)



# Velocity

df_peaks_summary <- df_peaks %>% 
  group_by(ppid, trial_num, feature) %>% 
  summarise(mean_opening_velocity = mean(velocity_thumb_finger[which(feature == "peak")]),
            mean_closing_velocity = mean(velocity_thumb_finger[which(feature == "trough")]))
ggplot(df_peaks_summary, aes(feature, mean_opening_velocity, color = ppid))+
  geom_boxplot()+
  geom_boxplot(aes(feature, mean_closing_velocity))+ 
  theme_bw()+
  labs(x='Peaks' ,y = "Velocity (mm/s)")  
#coord_cartesian(ylim = c(0, 1.5))




ggplot(df_peaks %>% dplyr :: filter(trial_num == 'Ah02T'), aes(x = time, y = velocity_thumb_finger, color = feature)) + 
  geom_point()  + 
  scale_x_continuous(limits = c(6, 9)) +
  #facet_wrap(~ trial_num)+ 
  theme_bw()+
  labs(x="Time (s)",y = "Velocity (mm/s)")  
#coord_cartesian(ylim = c(0, 1.5))

ggplot(df_peaks, aes(x = feature, y = 2*half_time_period, color = feature)) + 
  geom_boxplot() + 
  stat_summary(aes(group = trial_num), fun.data = mean_se, geom = "point", size = 2)+
  scale_y_continuous(limits = c(0, 0.5))+ 
  theme_bw()+
  labs(x='Peaks',y = "Time Period (s)")   #+
#facet_wrap(~ trial_num)


ggplot(df_peaks, aes(x = feature, y = freq, color = feature)) + 
  geom_boxplot() + 
  #stat_summary(aes(group = trial_num), fun.data = mean_se, geom = "point", size = 2)+
  scale_y_continuous(limits = c(0, 0.5))+ 
  theme_bw()+
  labs(x='Peaks',y = "Frequency (Hz")   #+
#facet_wrap(~ trial_num)

ggplot(df_peaks, aes(x = time, y = freq, group = ppid, color = feature)) + 
  geom_boxplot() + 
  #stat_summary(aes(group = trial_num), fun.data = mean_se, geom = "point", size = 2)+
  #scale_x_continuous(limits = c(6, 7))+ 
  #scale_y_continuous(limits = c(0, 0.2))+ 
  theme_bw()+
  labs(x='Time (s)',y = "Frequency (HZ)")   #+
#facet_wrap(~ trial_num)



ggplot(df_peaks%>% dplyr :: filter(trial_num == 'Ah02T'), aes(x = half_time_period, y = distance_thumb_finger)) + 
  geom_point(aes(color = feature)) + 
  stat_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, 0.5)) +
  facet_wrap(~ trial_num)

ggplot(df_peaks, aes(x = abs(velocity_thumb_finger), y = distance_thumb_finger, color = feature)) + 
  geom_point() + 
  facet_wrap(~ trial_num)


test_df <- df_peaks %>%
  dplyr::filter(time > 5, time < 10, trial_num == "Ah00T")


peaks <- test_df %>%
  pull(distance_thumb_finger) %>% 
  findpeaks()

test_df %>% 
  mutate(peak = )

df_movement_raw %>% 
  mutate(y_delta = distance_thumb_finger,
         y_diff_v = c(NA, diff(y_delta))) %>% 
  dplyr::filter(trial_num == 'Ah02T') %>% 
  ggplot(aes(x = time, y = y_diff_v)) + 
  geom_line() + 
  geom_line(aes(y = y_delta), color = 'red') + 
  #facet_wrap(~ trial_num, scales = "free") +
  scale_x_continuous(limits = c(6, 9))+ 
  theme_bw()+
  labs(x="Time (s)",y = "    Velocity    ,      Distance") 


#coord_cartesian(ylim = c(0, 1.5))



ggplot(df_movement_raw %>% dplyr::filter(trial_num =='Ah02T'), aes(time, distance_thumb_finger)) + 
  geom_path() + 
  scale_x_continuous(limits = c(2, 15))+ 
  theme_bw()+
  labs(x="Time (s)",y = "Distance_Thumb_Finger (mm)") + 
  coord_cartesian(ylim = c(0, 160))


ggplot(df_movement_raw, aes(time, distance_thumb_finger, color = interaction(ppid, trial_num))) + 
  geom_path() + 
  scale_x_continuous(limits = c(10, 12))


for(tn in unique(df_movement_raw$trial_num)){
  df_subset <- df_movement_raw %>% 
    dplyr::filter(trial_num == tn)
  
  #ttl <- sprintf("T: %s\nPPID: %s\nDIST: %scm", tn, df_subset$ppid[1], df_subset$distance_cm[1])
  
  (move_graph <- ggplot(df_subset,
                        aes(x = time, y = distance_thumb_finger, color = ppid)) + 
      geom_path(size = 1.2) + 
      scale_x_continuous(limits=c(10,11)) #+
    #coord_cartesian(ylim = c(0, 200))
  )
  fname <- sprintf('T%s.png',tn)
  ggsave(paste("distance_between",fname, sep = "/") , move_graph, type = "cairo-png")
}


for(tn in unique(df_movement_raw$trial_num)){
  df_subset <- df_movement_raw %>% 
    dplyr::filter(trial_num == tn)
  
  #ttl <- sprintf("T: %s\nPPID: %s\nDIST: %scm", tn, df_subset$ppid[1], df_subset$distance_cm[1])
  
  (move_graph <- ggplot(df_subset,
                        aes(x = time, y = finger_speed, color = ppid)) + 
      geom_path(size = 1.2)+ 
      scale_x_continuous(limits=c(10,11)) #+ 
    #coord_cartesian(ylim = c(0, 200))
  )
  fname <- sprintf('T%s.png',tn)
  ggsave(paste("finger_speed",fname, sep = "/") , move_graph, type = "cairo-png")
}

for(tn in unique(df_movement_raw$trial_num)){
  df_subset <- df_movement_raw %>% 
    dplyr::filter(trial_num == tn)
  
  #ttl <- sprintf("T: %s\nPPID: %s\nDIST: %scm", tn, df_subset$ppid[1], df_subset$distance_cm[1])
  
  (move_graph <- ggplot(df_subset,
                        aes(x = time, y = thumb_speed, color = ppid)) + 
      geom_path(size = 1.2)+ 
      scale_x_continuous(limits=c(10,11)) #+ 
    #coord_cartesian(ylim = c(0, 200))
  )
  fname <- sprintf('T%s.png',tn)
  ggsave(paste("thumb_speed",fname, sep = "/") , move_graph, type = "cairo-png")
}

for(tn in unique(df_movement_raw$trial_num)){
  df_subset <- df_movement_raw %>% 
    dplyr::filter(trial_num == tn)
  
  #ttl <- sprintf("T: %s\nPPID: %s\nDIST: %scm", tn, df_subset$ppid[1], df_subset$distance_cm[1])
  
  (move_graph <- ggplot(df_subset,
                        aes(x = time, y = finger_speed-thumb_speed, color = ppid)) + 
      geom_path(size = 1.2)+ 
      scale_x_continuous(limits=c(10,11)) #+ 
    #coord_cartesian(ylim = c(0, 200))
  )
  fname <- sprintf('T%s.png',tn)
  ggsave(paste("finger_speed_thumb_speed",fname, sep = "/") , move_graph, type = "cairo-png")
}
