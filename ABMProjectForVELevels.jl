module ABMProjectForVELevels


##SEEDING4 - This code allows us to set the epidemic to start on day 31 with
#500 seeding infections not amoung the vaccinated group, 100-0% vaccine protection
#that wanes over 60 days




## -- Name the Files You're Going to Output
function namedfile1(df,namefile)
  work_dir="//Users//nikas//Desktop//Batch5-30Day"
  name = work_dir*namefile
  CSV.write(name,df)
end
#work_dir would need to be changed based on Users, and folder pathway
#Generally names should be descriptive about the scenarios ran

function runABM()
##-- Make Name List Long Enough for the Number of Simulations (threads)
sim_number=1 #How many simulations to run
name = collect(1:sim_number)
name2 = collect((1+sim_number):(2*sim_number))
starttime = time()
#Random.seed!(383)


Threads.@threads for g in 1:sim_number

##-- Preset all variables and allocate space for vectors
D = 180  #Number of days in the season
N = 100000 # Population size
max_age = 1500 #currently unused
vacc_days= collect(1:30) #The days people get vaccinated CHANGE IF YOU WANT 60
vacc_freq=0.4  #The total precent of people who get vaccinated
init_freq=500.0/100000. #The total people INITIALLY infected (CDC says 2.4% baseline)
intr_day =1 #Day disease is introduced

beta_c = 5 # number of contacts
beta_I = 0.072 # infection chance 0.072
beta1 = beta_c*beta_I #this makes the standard beta for R0
gamma_1 = (1.0/5.0) #Recover in 5 days- this can be made into a distribution


vaccinated = zeros(Int64,N,3) #Save information about vaccinated (col 1 vac/not, col 2 when, col 3 protection if all/nothing)
infected = zeros(Int64,N)  # 1 infected, 0 not
recovered = zeros(Int64,N) # 1 recovered, 0 not
age =  rand(0:max_age,N)   #not currently in use, in case we want to incorporate death
total_daily_infection = zeros(Int64,D) #saves total new infections per day
exposures = zeros(Int64,D,6) # saves new exposures per day (v2v, v2u, u2v, u2u, inf, rec)
infections = zeros(Int64,D,4) # saves new infections per day (v2v, v2u, u2v, u2u)
infected_time = zeros(Int64,N) #when people get infected
recovered_time = zeros(Int64,N) #when people recover
ever_infected  = zeros(Int64,N) #if they were ever previously infected
original_infectors  = zeros(Int64,N) #if the were the intial seeders
ongoing = zeros(Int64,D) #current infection level
unvac_inf = zeros(Int64,D) #number of unvaccinated infected
vac_inf = zeros(Int64,D)#number of vaccinated infected
rec = zeros(Int64,D)#for saving recovery data
additional_vac =0
#-#############################################################################################################-#
# ------------------------------- RUN DAILY: AGING, VACCINATION, INFECTION ----------------------------------- #
#-#############################################################################################################-#

    for i in 1:D
      if (i == 1) #This is not being currently used but is for ageing population
          else

          age .+= 1
    end # for if it is a first day

##-- Run Seeding Infection
    if i== intr_day   # if it is the day the virus is introduced
      unvaccinated = collect(1:N)[(vaccinated[:,1] .==0) .& (ever_infected .== 0)] #find unvaccinated, uninfected people
      infect = sample(unvaccinated,round(Int,init_freq*N),replace = false) #randomly infect some people from the previous line
      infected[infect] .= 1 #they now have virions .= means scalar is passed to every element of array
      infected_time[infect] .= intr_day #save the day they got infected
      ever_infected[infect] .= 1 #marked as having gotten infected (even if they recover)
      original_infectors[infect] .= 1 #mark as the orginal seeders of infection

    end # for if int_day
##--Chose Whom to Vaccinate
      if i in vacc_days

        unvaccinated = collect(1:N)[(vaccinated[:,1] .==0) .& (ever_infected .== 0)] #find who is unvaccinated and not yet infected
         if length(unvaccinated) >= round(Int,vacc_freq/length(vacc_days)*N) #because of the rounding we typically get 39.99% vaccinated
          vaccinate = sample(unvaccinated,round(Int,vacc_freq/length(vacc_days)*N),replace = false)
          vaccinated[vaccinate,1] .= 1#mark as vaccinated
          vaccinated[vaccinate,2] .= i #save day of vaccination
          vaccinated[vaccinate,3] .= 1# protected (1) or not (0)(this changes if all or nothing)
          additional_vac_add =(vacc_freq/length(vacc_days)*N)-round(Int,vacc_freq/length(vacc_days)*N)
          additional_vac=additional_vac+additional_vac_add
          if additional_vac >1
            vaccinate = sample(unvaccinated,round(Int,1),replace = false)
            vaccinated[vaccinate,1] .= 1#mark as vaccinated
            vaccinated[vaccinate,2] .= i #save day of vaccination
            vaccinated[vaccinate,3] .= 1
            additional_vac=0
          end
          #NEXT 2 Lines: All or Nothing 80% (20% of vaccinated would have no protection)
          #unprotect1 =sample(vaccinate,round(0.2*length(vaccinate)),replace = false)
          #vaccinated[unprotect1,3].=0

          #Optional Test Case: To ensure correct, set protection to 0, Cox VE should give 0% VE
        end # of length of unvaccinated is short enough
        end # of if i in vacc days
      ###########



#-#############################################################################################################-#
# ------------------------------------------ INFECT PEOPLE ZONE ----------------------------------------------- #
#-#############################################################################################################-#

##--Protection Function For Vaccinated


leaking = collect(0:0.01666:1) #100 to 0 protection over 60 days
#leaking=ones(180)*0.2 #80% Proection
    #  leaking = collect(0.2:0.005:0.5) #80 to 50 protection over 365 days
      ##Optional Test Case: Make leaking all ones (no protection) should recieve all 0 Cox
      #leaking =ones(60)
infect_probsv2v = append!(leaking, ones(121)) #append so it is 80% for the remaining time (this line is more important for waning cases)
      infect_probsu2v = infect_probsv2v  # infection by unvaccinated and vaccinated assumed to be the same
      ##-- Protection Function for Unvaccinated
      infect_probsv2u = 1
      infect_probsu2u = 1
      trans_prob_fcn = beta_I #This incorporates beta in the chance to get infected
      gen_contacts =  Array{Int64}(undef,beta_c,1) #generate 5 contacts

##--Save data for contacts and exposures everyday (initalize at 0)
      vac2vac_inf = 0
      vac2unv_inf = 0
      unv2vac_inf=0
      unv2unv_inf=0
      vac2vac_exp = 0
      vac2unv_exp=0
      unv2vac_exp=0
      unv2unv_exp=0
      inf_exp=0
      rec_exp=0
##--Checking all people
            for j in 1:N

        #IF YOU ARE INFECTED
              if infected[j]==1
                    ever_infected[j] = 1 #ensuring they are correctly marked (this is a failsafe)
           #if infected for 5 days, recover, previously 5 =chance[j] so that it was variable
                if infected_time[j] == (i-5)
                    infected[j]=0 #now recovered after 5 days
                    recovered[j]=1#now recovered after 5 days
                    recovered_time[j]=i #save when they got recovered
                  else #otherwise if still infectious, generate contacts

                   gen_contacts[:,1] = sample(1:N, beta_c, replace = false) #generates 5 random contacts
                   for index in gen_contacts[:,1]
                    acquire = rand() #generate random acquire number (if less then get infected, later lines)
                    if recovered[index] == 1 || infected[index] == 1 #if recovered or currently infected, contact cannot get reinfected
                      if recovered[index] == 1
                        rec_exp+=1
                      end
                      if infected[index] == 1
                        inf_exp+=1
                      end
                    else
                      if vaccinated[j,1]==1 && vaccinated[index,1]==1 && vaccinated[index,3] ==1 #if both people vaccinated

                      if  acquire <= infect_probsv2v[i-vaccinated[index,2]+1]*trans_prob_fcn #chance to acquire based on when vaccinated for potential infectee
                        infected[index] = 1 #if acquired label as infected
                        ever_infected[index] = 1 #label as ever infected (for after they recover)
                        infected_time[index] = i #give time of infection
                         vac2vac_inf += 1 #add one to the infected that day

                       else #if they didnt get infected count as an exposure
                         vac2vac_exp += 1
                        end # of if with acquire condition
                      end #end vac 2 vac infection condition

                      ##-- This is just for all or nothing vaccinated to vaccinated (not currently in use)
                      if vaccinated[j,1]==1 && vaccinated[index,1]==1 && vaccinated[index,3] ==0

                      if  acquire <= infect_probsu2u*trans_prob_fcn
                        infected[index] = 1
                        ever_infected[index] = 1
                        infected_time[index] = i
                         vac2vac_inf += 1

                         else

                         vac2vac_exp += 1
                        end # of if with acquire condition
                      end #end all or nothing vac to vac
##--Vaccinated people infecting unvaccinated people
                     if vaccinated[j,1]==1 && vaccinated[index,1]==0

                      if  acquire <= infect_probsv2u*trans_prob_fcn
                        infected[index] = 1
                        ever_infected[index] = 1
                        infected_time[index] = i
                         vac2unv_inf += 1

                      else #if uninfected count exposure
                         vac2unv_exp +=1
                    end # of if with aqcuire
                  end #end vac 2 unv infection
##-- Unvaccinated People Infecting Vaccinated People (who have protection)
                  if vaccinated[j,1]==0 && vaccinated[index,1]==1 &&vaccinated[index,3]==1
                      if  acquire <= infect_probsu2v[i-vaccinated[index,2]+1]*trans_prob_fcn
                        infected[index] = 1
                        ever_infected[index] = 1
                        infected_time[index] = i
                         unv2vac_inf +=1
                     else
                     #if uninfected count the exposure
                        unv2vac_exp += 1
                    end
                  end
##-- For all-or-nothing unvaccinated infecting vaccinated (not currently in use)
                  if vaccinated[j,1]==0 && vaccinated[index,1]==1 &&vaccinated[index,3]==0
                     #Infector Un, Infecctee Vac
                      #if(length(which(acquire<=infect.probsu2v[i-vaccinated[j,2]]*trans.prob.fcn))!=0){
                      if  acquire <= infect_probsu2u*trans_prob_fcn
                        infected[index] = 1
                        ever_infected[index] = 1
                        infected_time[index] = i
                         unv2vac_inf +=1

                     else
                     #if uninfected count the exposure
                        unv2vac_exp += 1
                      #}
                    end
                  end
##--Unvaccinated to unvaccinated Infections
                  if vaccinated[j,1]==0 && vaccinated[index,1]==0
                      if  acquire <= infect_probsu2u*trans_prob_fcn
                        infected[index] = 1
                        ever_infected[index] = 1
                        infected_time[index] = i
                         unv2unv_inf += 1
                        else
                          #if uninfected count the exposure
                          unv2unv_exp += 1
                      #}
                    end # of aquire condition
                  end # of unvac unvac combo
                  #}
                end # if recovered or vaccinated condition
                #}
              end # of cycle through potential infected indexes
              #}
            end # of if with time condition
           #}
         end # of if with is infected condition
        #}
      end # cycle through all people
  #-- END INFECTING/EXPOSING PEOPLE --#

## ------------------- SAVE INFECTED/EXPOSED DAILY INFORMATION IN AN EASY TO USE FORMAT ---------------------------- ##
      #first person is always infected, second person is who they attempt to pass it to
      #this is needed for Halloran VE calculations
      exposures[i,1] = vac2vac_exp
      exposures[i,2] = vac2unv_exp #these are infected vac exposing to unvac uninfected people
      exposures[i,3] = unv2vac_exp#these are infected unvac exposing to vac uninfected people
      exposures[i,4] = unv2unv_exp
      exposures[i,5] = inf_exp
      exposures[i,6] = rec_exp
      infections[i,1] = vac2vac_inf
      infections[i,2] = vac2unv_inf
      infections[i,3]= unv2vac_inf
      infections[i,4] = unv2unv_inf
       total_daily_infection[i] = vac2vac_inf+vac2unv_inf+unv2vac_inf+unv2unv_inf
       ongoing[i] = sum(infected)

##-- For each fday count up total of vaccinated infections, unvaccinated infections, and recoveries
        for k in 1:N
          if  infected[k]==1
              if vaccinated[k,1] == 0
                unvac_inf[i] = unvac_inf[i]+1
              else
                vac_inf[i] = vac_inf[i]+1
              end #condition for vaccinated
        end #condition for infected
          if recovered[k] == 1
             rec[i] = rec[i]+1
          end

       end # of k cycle
##-- Daily Cycle of Infection/Exposure Ends here
  end # where daily cycle ends


time_points = collect(1:D)
##-- Make a test for Cox model where it should give the exact opposite result
     vac_status = vaccinated[:,1] #ever vaccinated, not necessarily actively vaccinated

    vac_status2 = ones(Int64,N) # array of ones


    #for(h in 1:N){
    for h in 1:N
      #if(vac.status[h]==1){
       if vac_status[h] ==1
         #vac.status2[h]=0
         vac_status2[h] = 0
      #}
    end # of if statement
    # }
  end # for h cycle
  #--END THE COX TEST ZONE
##--CENSOR DATA FOR COX MODEL
    vac_time = vaccinated[:,2]

    id = collect(1:N)
    for n in 1:N
      #if(infected_time[n]== 0){
      if infected_time[n]== 0
#Censor Data for those who never got infected (having it be 0 would be misinterpreted in COX)
        infected_time[n]= D # Change to D?
      end # of if
     #}
   end # of for loop through n

##--CREATE DATA FRAMES TO SAVE OUT
    cum_AR = sum(ever_infected)
    percent = cum_AR/N
    mydf = DataFrame(id=id,vac_status = vac_status, ever_infected = ever_infected, infected_time = infected_time,
                     vac_time = vac_time, vac_status2 = vac_status2, original_infectors=original_infectors) #For Cox Model
    mydf2 = DataFrame(time_points=time_points,vac_inf = vac_inf,unvac_inf = unvac_inf, vac2vac_exp = exposures[:,1],
    vac2unv_exp = exposures[:,2], unv2vac_exp=exposures[:,3],unv2unv_exp=exposures[:,4],vac2vac_inf=infections[:,1],
    vac2unv_inf=infections[:,2],unv2vac_inf=infections[:,3],unv2unv_inf=infections[:,4], inf_exp=exposures[:,5], rec_exp=exposures[:,6])#For All Other VE
##--Save Data Out
    namedfile1(mydf,"$(g).csv")
    namedfile1(mydf2,"$(g+sim_number).csv") #Note that the +# should be equal to however many threads

end #End the simulations


timing2=time()-starttime #calculate how long it took (generally <3min)

end

end
