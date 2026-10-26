######################
# model HDFTS of CDFs
######################

# load R package

require(HMDHFDplus)
require(demography)
require(ftsa)
require(hdftsa)

# read Japanese Subnational Human Mortality Data

state = c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima","Ibaraki", "Tochigi", 
          "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", 
          "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", 
          "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", "Tokushima", 
          "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", 
          "Kagoshima", "Okinawa")

# change file directory

ages = 0:110
n_age = length(ages)
years = 1947:2023
n_year = length(years)

# for Japanese data, let us consider years from 1975 to 2023, ages from 0 to 94 in single age group, last age 95+ 

chosen_ages = 0:95
chosen_years = 1975:2023
years_index = 29:77

# raw and smooth demogdata

female_prefecture_rate = male_prefecture_rate = 
smooth_female_prefecture_rate = smooth_male_prefecture_rate = list()
prefecture_n_year = vector("numeric", length(state))
for(ik in 1:length(state))
{
    if(ik < 10)
    {
        code = paste("0",ik,sep="")
    }
    else
    {
        code = as.character(ik)
    }
    
    # females
    
    data_rate_F = matrix(readJMDweb(code, item = "Mx_1x1", fixup = FALSE)[,3], n_age, )
    data_expo_F = matrix(readJMDweb(code, item = "Exposures_1x1", fixup = FALSE)[,3], n_age, )
    
    demo_F = demogdata(data = data_rate_F, pop = data_expo_F, ages = ages, years = tail(years, ncol(data_rate_F)), 
                       type = "mortality", label = state[ik], name = "F")
    female_prefecture_rate[[ik]] = extract.ages(demo_F, chosen_ages)
    smooth_female_prefecture_rate[[ik]] = extract.ages(smooth.demogdata(demo_F), chosen_ages)
    
    # males
    
    data_rate_M = matrix(readJMDweb(code, item = "Mx_1x1", fixup = FALSE)[,4], n_age, )
    data_expo_M = matrix(readJMDweb(code, item = "Exposures_1x1", fixup = FALSE)[,4], n_age, )
    
    demo_M = demogdata(data = data_rate_M, pop = data_expo_M, ages = ages, years = tail(years, ncol(data_rate_M)), 
                       type = "mortality", label = state[ik], name = "M")
    male_prefecture_rate[[ik]] = extract.ages(demo_M, chosen_ages)
    smooth_male_prefecture_rate[[ik]] = extract.ages(smooth.demogdata(demo_M), chosen_ages)
    
    prefecture_n_year[ik] = ncol(data_rate_F)
    print(ik); rm(ik); rm(code); rm(data_rate_F); rm(data_expo_F); rm(data_rate_M); rm(data_expo_M)
    rm(demo_F); rm(demo_M)
}

# smoothing noise

smooth_female_prefecture_rate_noise = smooth_male_prefecture_rate_noise = array(NA, dim = c(96, 49, 47))
for(ik in 1:47)
{
    if(ncol(smooth_female_prefecture_rate[[ik]]$obs.var$F) == 77)
    {
        smooth_female_prefecture_rate_noise[,,ik] = smooth_female_prefecture_rate[[ik]]$obs.var$F[,29:77]
        smooth_male_prefecture_rate_noise[,,ik]   = smooth_male_prefecture_rate[[ik]]$obs.var$M[,29:77]
    }
    else if(ncol(smooth_female_prefecture_rate[[ik]]$obs.var$F) == 51)
    {
        smooth_female_prefecture_rate_noise[,,ik] = smooth_female_prefecture_rate[[ik]]$obs.var$F[,3:51]
        smooth_male_prefecture_rate_noise[,,ik]   = smooth_male_prefecture_rate[[ik]]$obs.var$M[,3:51]
    }
    rm(ik)
}

# plot Okinawa

# savefig("Fig_1a", width = 12, height = 10, toplines = 0.8, type = "png")
# plot(female_prefecture_rate[[47]], main = "Okinawa female data")
# dev.off()
# 
# savefig("Fig_1b", width = 12, height = 10, toplines = 0.8, type = "png")
# plot(male_prefecture_rate[[47]], main = "Okinawa male data")
# dev.off()

# convert list to matrix

female_prefecture_rate_array = male_prefecture_rate_array = 
female_prefecture_rate_array_smooth = male_prefecture_rate_array_smooth = array(NA, dim = c(length(chosen_ages), 47, length(years_index)))

female_prefecture_rate_combo = male_prefecture_rate_combo = 
female_prefecture_rate_combo_smooth = male_prefecture_rate_combo_smooth = matrix(NA, length(chosen_ages), (length(years_index) * 47))

female_prefecture_rate_combo = male_prefecture_rate_combo = 
female_prefecture_rate_combo_smooth = male_prefecture_rate_combo_smooth = NULL
for(ik in 1:length(state))
{
    if(ncol((female_prefecture_rate[[ik]]$rate$F)) == 77)
    {
        # raw
      
        female_prefecture_rate_array[,ik,] = (female_prefecture_rate[[ik]]$rate$F)[,years_index]
        male_prefecture_rate_array[,ik,]   = (male_prefecture_rate[[ik]]$rate$M)[,years_index]
        
        female_prefecture_rate_combo = cbind(female_prefecture_rate_combo, female_prefecture_rate_array[,ik,])
        male_prefecture_rate_combo = cbind(male_prefecture_rate_combo, male_prefecture_rate_array[,ik,])
        
        # smooth
        
        female_prefecture_rate_array_smooth[,ik,] = (smooth_female_prefecture_rate[[ik]]$rate$F)[,years_index]
        male_prefecture_rate_array_smooth[,ik,]   = (smooth_male_prefecture_rate[[ik]]$rate$M)[,years_index]
        
        female_prefecture_rate_combo_smooth = cbind(female_prefecture_rate_combo_smooth, female_prefecture_rate_array_smooth[,ik,])
        male_prefecture_rate_combo_smooth = cbind(male_prefecture_rate_combo_smooth, male_prefecture_rate_array_smooth[,ik,])
    }
    else if(ncol((female_prefecture_rate[[ik]]$rate$F)) == 51)
    {
        # raw
      
        female_prefecture_rate_array[,ik,] = (female_prefecture_rate[[ik]]$rate$F[,3:51])
        male_prefecture_rate_array[,ik,] = (male_prefecture_rate[[ik]]$rate$M[,3:51])
        
        female_prefecture_rate_combo = cbind(female_prefecture_rate_combo, female_prefecture_rate_array[,ik,])
        male_prefecture_rate_combo = cbind(male_prefecture_rate_combo, male_prefecture_rate_array[,ik,])
        
        # smooth
        
        female_prefecture_rate_array_smooth[,ik,] = (smooth_female_prefecture_rate[[ik]]$rate$F[,3:51])
        male_prefecture_rate_array_smooth[,ik,] = (smooth_male_prefecture_rate[[ik]]$rate$M[,3:51])
        
        female_prefecture_rate_combo_smooth = cbind(female_prefecture_rate_combo_smooth, female_prefecture_rate_array_smooth[,ik,])
        male_prefecture_rate_combo_smooth = cbind(male_prefecture_rate_combo_smooth, male_prefecture_rate_array_smooth[,ik,])
    }
    print(ik); rm(ik)
}

