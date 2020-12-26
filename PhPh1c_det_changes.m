%PhPh1c_det_changes function file. It is the second file accessed during 
%the simulation and the first file accessed by the PhPh1c_main file. 
%Requires two row vectors as input (The function will prompt the user to 
%input these row vectors). First row vector is the det_time vector. For 
%this vector, the user must input a row vector consisting of the times at 
%which deterministic changes must be made to thevsystem. 
%The second row vector is the det_changes vector. The user must input a row 
%vector consisting of the changes that need to be made at the
%times entered in the det_time vector.  
%Note 1: Both vectors must have the same number of elements, else the 
%larger vector will be truncated.
%Note 2: The times entered need not be in any order, the code will sort the
%vector.
%Note 3: If no changes are to be made to the queue at any time during the
%simulation, enter empty matrix ([]) when prompted in both cases
%Note 4: If the user does not wish to start the system in empty and idle
%state (default setting), then the number of entities in the queue at the
%start must be entered (Enter 0 as one of the elements in the det_time 
%vector and the corresponding values in the det_change vector will be taken
%as the initial number of entities in queue).                              
%Note 5: There is no provision to ensure that only non negative values are
%entered as the initial number of entities in the queue. Hence the user
%must take care not to enter a negative value for the initial number of
%entities in the queue. Similarly, there is no provision to ensure 
%negative values for time entered in det_time vector will be removed. 
%The user must ensure this. However, negative values are allowed as 
%deterministic changes, albeit for times that are not equal to the
%start_time.
%Note 6: At times other than start time, there is no restriction for the
%value of the change to be made with respect to the number of entities in
%the queue. But any value greater than c or lesser than -c will be taken as
%c and -c respectively. 
%Note 7: The array z that is returned consists of the following 7 main 
%vectors (These will be referred to as data vectors henceforth):
%t : simulation time values
%E0_feq : 0th PMDE values calculated using the Kolmogorov Forward
%equations (KFE)
%E1_feq : 1st moment values calculated using the KFE
%E2_feq : 2nd moment values calculated using the KFE
%E0_pmde : 0th PMDE values calculated using the PE approximation
%E1_pmde: 1st PMDE values calculated using the PE approximation
%E2_pmde: 2nd PMDE values calculated using the PE approximation
%PE approximation stands for Polya Eggenberger approximation. 
function [z]=PhPh1c_det_changes()

global c max_arrival_phases max_service_phases det_changes det_time;
global start_time end_time;

%%%%%%%%%%%%%%%%%%% Input det_time and det_changes %%%%%%%%%%%%%%%%%%%%%%%%
%det_time denotes deterministic change time
%det_changes denotes deterministic changes

fprintf('\nIf no change is entered for time 0, then the system will start as empty and idle\n');

%Prompt for row vector that will specify the time at which changes need to
%be made.
det_time=input('\nEnter the times at which changes must be made as a row vector: \n');
%Prompt for row vector that will specify the changes corresponding to the
%times given the det_time vector that need to be made. 
det_changes=input('\nEnter the changes corresponding to the time entered above as a row vector: \n');

%Both row vectors are converted to column vectors purely for
%computational purposes. The user is asked to enter these two vectors as
%row vectors due to convenience of entering data in a row. 
det_time=det_time';
det_changes=det_changes';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If both vectors don't have the same number of elements, truncate the
% larger vector. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_det_time=size(det_time);
s_det_changes=size(det_changes);

if s_det_time<s_det_changes
    det_changes=det_changes(1:s_det_time,1);
end

if s_det_time>s_det_changes
    det_time=det_time(1:s_det_changes,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% Sort arrays and remove incorrect values %%%%%%%%%%%%%%%%%%%%%%

%The code appends the start_time and end_time to the det_time matrix if no
%changes are made to the queue at these times. This is purely for
%computational purposes

%This section of the code will sort the vectors and remove incorrect values,
%i.e. values in the det_time array that are greater than end_time or lesser
%than start_time. If there is a det_time element that is incorrect, then
%its corresponding det_changes element will be removed as well. 

%Sort both vectors only if they are not empty. 
if sum(size(det_time))>0
    
[det_time,I]=sort(det_time); 

count=size(det_time);
count=count(1,1);
dummy=det_changes;

for i=1:count
    index=I(i,1);
    det_changes(i,1)=dummy(index,1);
end

clear dummy;

%Store sorted det_time and det_changes vectors in dummy vectors, clear
%det_time and det_change. Then insert values from the dummy vectors to
%det_time and det_change vectors if they are not incorrect values. 

dummy1=det_time;     
dummy2=det_changes;

max_index=count;   %count = original length of det_time vector.

det_time=0;
det_changes=0;
count=1;

for i=1:max_index
    if dummy1(i,1)>=start_time 
        if dummy1(i,1)<=end_time
            det_time(count,1)=dummy1(i,1);
            det_changes(count,1)=dummy2(i,1);
            count=count+1;
        end
    end
end

count=count-1;

dummy1=det_time;
dummy2=det_changes;

%check if the user has input start_time as one of the values. If so, do not
%append start_time to det_time vector. Else append start_time to det_time
%vector. 
if dummy1(1,1)==start_time
    det_time(1:count,1)=dummy1(1:count,1);
    det_changes(1:count,1)=dummy2(1:count,1);
else
    det_time(1,1)=start_time;
    det_changes(1,1)=0;
    det_time(2:count+1,1)=dummy1(1:count,1);
    det_changes(2:count+1,1)=dummy2(1:count,1);
end

%check if the user has input end_time as one of the values. If so, do not
%append end_time to det_time vector. Else append end_time to det_time
%vector. 

dummy=size(det_time);
dummy=dummy(1,1);

if dummy1(count,1)~=end_time
    det_time(dummy+1,1)=end_time;
    det_changes(dummy+1,1)=0;
end

clear dummy dummy1 dummy2;
    
else
    %If the user entered [] for both prompts, then det_time array consists
    %of start_time and end_time and det_changes array consists of 0 and 0. 
    det_time(1,1)=start_time;
    det_time(2,1)=end_time;
    det_changes(1,1)=0;
    det_changes(2,1)=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Det_changes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This part of the code will run the ode45 in a loop the required number of
%times. Each run goes from the time of last change (In the first run, this
%is the start_time) to the time of next change (If no changes are made in
%the queue at any time during the simulation then this is the end_time).Two 
%sets of values are calculated: First set calculates the 0th, 1st and 2nd 
%PMDE by using the exact values of state probabilities (Obtained by 
%evaluating the Kolmogorov Forward Equations)and the second set calculates 
%the 0th, 1st and 2nd PMDE by using the Polya Eggenberger distribution 
%to approximate the required probabilites at change times. For more
%information on the sets of equations, refer to PhPh1c_function file. 


%At the end of each run, values calculated will be concatenated with the
%corresponding main vector. For example: if a change is made at time = 2,
%where start_time = 0 and end_time =10. Then values are first calculated
%for the time period from 0 to 2. The values are then stored in their
%corresponding data vectors (Refer to the information at the start of this
%file to know what are the data vectors). The next run is from time 2 to 
%time 10. Since we have a new value at time 2 due to the change, the
%current values will be concatenated to the last element in each data
%vector, thus erasing the the last time step for the previous run. Note:
%ode45 ensures that the last time step is the end time specified. Thus the
%last value of time is almost always the last value specified in the ode45
%paramter tspan. 

length_det_time=size(det_time);
length_det_time=length_det_time(1,1);

opts=odeset('RelTol',1e-6);

length_t=1;


%Declare the data vectors. 
t=0;
E0_pmde=0;
E1_pmde=0;
E2_pmde=0;

E0_feq=0;
E1_feq=0;
E2_feq=0;

%Run the ode45 in a loop

for i=1:length_det_time-1 
    
    s_time=det_time(i,1);  %start time for the current run
    e_time=det_time(i+1,1); %end time for the current run
    duration=e_time-s_time;  %duration of the current run 
    tspan=linspace(s_time,e_time,duration*100); %tspan for the current run
    change=det_changes(i,1);
    
    %This section of the code relies heavily on the functioning of the
    %PhPh1c_initial_conditions file. It would be fair to say the logic
    %behind the deterministic changes to the Ph(t)/Ph(t)/1/c queue is in
    %that file. The user must refer to the function file
    %PhPh1c_initial_conditions file to understand the whole code in
    %general.
    
    if i==1
        ics=PhPh1c_initial_conditions(0,change,0);
        [time,x]=ode45(@PhPh1c_function,tspan,ics,opts);
        length_x=size(x);
        length_x=length_x(1,1);
        
        next_ics=(x(length_x,:))';
    else
        ics=PhPh1c_initial_conditions(1,change,next_ics);
        [time,x]=ode45(@PhPh1c_function,tspan,ics,opts);
        length_x=size(x);
        length_x=length_x(1,1);
        
        next_ics=(x(length_x,:))';
    end
    
    
    %At end of each run, the values calculated using the PhPh1c_function
    %file and ode45 will be concatenated to the required vectors. 
    
    
    E0_feq_dummy=0;
    E1_feq_dummy=0;
    E2_feq_dummy=0;
    
    %Refer to the PhPh1c_function file for more information on how pmde
    %values are calculated and what the array x means.  
    
    
    %we use dummy variables for each data vector, to make required 
    %calculations and then concatenate the dummy vectors to the data 
    %vectors. 
for n=1:c
    E0_feq_dummy=E0_feq_dummy+sum((x(:,max_arrival_phases+...
               (n-1)*max_arrival_phases*max_service_phases+1:...
               n*max_arrival_phases*max_service_phases))')';
    
    
    E1_feq_dummy=E1_feq_dummy+n*sum((x(:,max_arrival_phases+...
               (n-1)*max_arrival_phases*max_service_phases+1:...
               n*max_arrival_phases*max_service_phases))')';
           
    
    E2_feq_dummy=E2_feq_dummy+(n^2)*sum((x(:,max_arrival_phases+...
               (n-1)*max_arrival_phases*max_service_phases+1:...
               n*max_arrival_phases*max_service_phases))')';
    
end       

E0_pmde_dummy=(sum(x(:,max_arrival_phases+...
                   c*max_arrival_phases*max_service_phases+1:...
                   max_arrival_phases+...
                   (c+1)*max_arrival_phases*max_service_phases)'))';    

E1_pmde_dummy=(sum(x(:,max_arrival_phases+...
                   (c+1)*max_arrival_phases*max_service_phases+1:...
                   max_arrival_phases+...
                   (c+2)*max_arrival_phases*max_service_phases)'))';    


E2_pmde_dummy=(sum(x(:,max_arrival_phases+...
                   (c+2)*max_arrival_phases*max_service_phases+1:...
                   max_arrival_phases+...
                   (c+3)*max_arrival_phases*max_service_phases)'))';    

     
    %Concatenate dummy vectors to main vectors.            
    t(length_t:length_t+length_x-1,1)=time(1:length_x,1);
    E0_feq(length_t:length_t+length_x-1,1)=E0_feq_dummy(1:length_x,1);
    E1_feq(length_t:length_t+length_x-1,1)=E1_feq_dummy(1:length_x,1);
    E2_feq(length_t:length_t+length_x-1,1)=E2_feq_dummy(1:length_x,1);
    E0_pmde(length_t:length_t+length_x-1,1)=E0_pmde_dummy(1:length_x,1);
    E1_pmde(length_t:length_t+length_x-1,1)=E1_pmde_dummy(1:length_x,1);
    E2_pmde(length_t:length_t+length_x-1,1)=E2_pmde_dummy(1:length_x,1);
    length_t=size(t);
    length_t=length_t(1,1);
    
    clear length_x time x;
    clear E0_feq_dummy E1_feq_dummy E2_feq_dummy;
    clear E0_pmde_dummy E1_pmde_dummy E2_pmde_dummy;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The last run is done manually to check if any changes have been made at
%end_time. If so, we can simply use the PhPh1c_initial_conditions file to
%obtain the required values. 
if det_changes(i+1,1)~=0
    change=det_changes(i+1,1);
    x=PhPh1c_initial_conditions(1,change,next_ics);
    length_x=size(x);
    length_x=length_x(1,1);
        

    E0_feq_dummy=0;
    E1_feq_dummy=0;
    E2_feq_dummy=0;
    
    
for n=1:c
    E0_feq_dummy=E0_feq_dummy+sum((x(:,max_arrival_phases+...
               (n-1)*max_arrival_phases*max_service_phases+1:...
               n*max_arrival_phases*max_service_phases))')';
    
    
    E1_feq_dummy=E1_feq_dummy+n*sum((x(:,max_arrival_phases+...
               (n-1)*max_arrival_phases*max_service_phases+1:...
               n*max_arrival_phases*max_service_phases))')';
           
    
    E2_feq_dummy=E2_feq_dummy+(n^2)*sum((x(:,max_arrival_phases+...
               (n-1)*max_arrival_phases*max_service_phases+1:...
               n*max_arrival_phases*max_service_phases))')';
    
end       

E0_pmde_dummy=(sum(x(:,max_arrival_phases+...
                   c*max_arrival_phases*max_service_phases+1:...
                   max_arrival_phases+...
                   (c+1)*max_arrival_phases*max_service_phases)'))';    

E1_pmde_dummy=(sum(x(:,max_arrival_phases+...
                   (c+1)*max_arrival_phases*max_service_phases+1:...
                   max_arrival_phases+...
                   (c+2)*max_arrival_phases*max_service_phases)'))';    


E2_pmde_dummy=(sum(x(:,max_arrival_phases+...
                   (c+2)*max_arrival_phases*max_service_phases+1:...
                   max_arrival_phases+...
                   (c+3)*max_arrival_phases*max_service_phases)'))';    

    
    t(length_t:length_t+length_x-1,1)=time(1:length_x,1);
    E0_feq(length_t:length_t+length_x-1,1)=E0_feq_dummy(1:length_x,1);
    E1_feq(length_t:length_t+length_x-1,1)=E1_feq_dummy(1:length_x,1);
    E2_feq(length_t:length_t+length_x-1,1)=E2_feq_dummy(1:length_x,1);
    E0_pmde(length_t:length_t+length_x-1,1)=E0_pmde_dummy(1:length_x,1);
    E1_pmde(length_t:length_t+length_x-1,1)=E1_pmde_dummy(1:length_x,1);
    E2_pmde(length_t:length_t+length_x-1,1)=E2_pmde_dummy(1:length_x,1);
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Store data vectors in the array named z. 

z(:,1)=t;
z(:,2)=E0_feq;
z(:,3)=E1_feq;
z(:,4)=E2_feq;
z(:,5)=E0_pmde;
z(:,6)=E1_pmde;
z(:,7)=E2_pmde;

end