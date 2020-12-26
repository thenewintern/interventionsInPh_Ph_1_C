%%%%%%%%%%%%%%%%%%%%% Main File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the file that needs to be run on the command line. To run this
%file type the file name (PhPh1c_main) on the command line. The program
%will prompt the user twice for two inputs. For more information on the
%inputs required by the program, refer to the function file
%PhPh1c_det_changes.
clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Declare global variables. 
global c max_arrival_phases max_service_phases det_changes det_time;
c=10;
max_arrival_phases=10;
max_service_phases=10;
%Variable usage and information:
%c : the total capacity of the queue. Thus at any given time, at most c-1
%customers are waiting for service and any customers that arrive when 
%c-1 customers are waiting for service are cleared.
%max_arrival_phases : the maximum number of arrival phases that are allowed
%in this simulation.
%max_service_phases : the maximum number of service phases that are allowed
%in this simulation.
%This code is designed to accept a maximum of 10 arrival phases and 10
%service phases. This ceiling can be increased by altering values in the
%PhPh1c_qparm file. Note, when the number of phases in the arrival process
%or service process are less than 10, no modifications need to be made to
%any part of the code. If the user wishes to run the program for a number
%of arrival or service phases greater than 10, changes must be made in the
%file PhPh1c_qparm. Refer to the function file PhPh1c_qparm for more
%information on these changes. 
%det_changes : stores the number of changes made to entities in the queue 
%det_time : stores the times at which these changes need to be made
%det_changes and det_time are the two inputs required by the user (In the
%form of row vectors). For more information on the input, refer to
%PhPh1c_det_changes function file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global start_time end_time;
start_time=0;
end_time=100;
%start_time denotes the starting time value for the simulation and
%similarly end_time denotes the end time of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[z]=PhPh1c_det_changes();
%Refer to the PhPh1c_det_changes function file for more information on this
%function and the array returned (z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z is an array of 7 column vectors (Which will be referred to as data
%vectors in the file PhPh1c_det_changes).
%Vector name and usage:
%t : stores the values of time at which the ode45 calculates values of the
%differential equations (Time steps)
%E0_feq : 0th PMDE values, calculated at each time step using the
%Kolmogorov Forward equations
%E1_feq : 1st PMDE values, calculated at each time step using the
%Kolmogorov forward equations
%E2_feq : 2nd PMDE values, calculated at each time step using the
%Kolmogorv forward equations
%E0_pmde : 0th PMDE values, calculated at each time step using the
%Partial moment differential equations (PMDEs) and Polya Eggenbeger (PE) 
%approximations.                                                           
%E1_pmde : 1st PMDE values, calculated at each time step using the PMDEs
%and PE approximations
%E2_pmde : 2nd PMDE values, calculated at each time step using the PMDEs
%and PE approximations. 

t=z(:,1);
E0_feq=z(:,2);
E1_feq=z(:,3);
E2_feq=z(:,4);
E0_pmde=z(:,5);
E1_pmde=z(:,6);
E2_pmde=z(:,7);
Var_feq=E2_feq-(E1_feq.*E1_feq);
Var_pmde=E2_pmde-(E1_pmde.*E1_pmde);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Once the data is obtained, calculate the error between the actual values
%and the approximations. 
            
E0_error=E0_feq-E0_pmde; %Difference between actual E0 values and approximate E0 values              
E1_error=E1_feq-E1_pmde; %Difference between actual E1 values and approximate E1 values              
Var_error=Var_feq-Var_pmde; %Difference between actual E2 values and approximate E2 values              
E2_error=E2_feq-E2_pmde;
epsilon=1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following section calculates Max absolute error and the corresponding
%relative error along with Max relative error and the corresponding
%absolute error for the first and second moments. Note: This data is not
%automatically presented by the code in the output. The variables will be
%present in the workspace and must be accessed either by directly checking
%its value (Given in the work space) or typing its name. 

%To understand the following section, note that there may be multiple
%values that have the same max absolute error which will ensure we have
%multiple corresponding relative errors. In order to counter this, we will
%consider only the maximum corresponding relative error. The logic remains
%the same for max relative error. 

%First Moment
dummy=abs(E1_error);

Max=max(dummy)-epsilon;
dummy=floor(dummy/Max);

index=find(dummy);

max_absolute_error_E1=max(abs(E1_error(index,1)));
max_absolute_error_corresponding_relative_error_E1=max(...
                          (abs(E1_error(index,1))./E1_feq(index,1))*100);

                               
dummy=abs(Var_error);

Max=max(dummy)-epsilon;
dummy=floor(dummy/Max);

index=find(dummy);

max_absolute_error_Var=max(abs(Var_error(index,1)));
max_absolute_error_corresponding_relative_error_Var=max(...
                        (abs(Var_error(index,1))./Var_feq(index,1))*100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                               

%Second Moment
dummy=(abs(E1_error)./E1_feq)*100;

Max=max(dummy)-epsilon;
dummy=floor(dummy/Max);

index=find(dummy);

max_relative_error_E1=max((abs(E1_error(index,1))./E1_feq(index,1)))*100;
max_relative_error_corresponding_absolute_error_E1=...
                                               max(abs(E1_error(index,1)));


dummy=(abs(Var_error)./Var_feq)*100;

Max=max(dummy)-epsilon;
dummy=floor(dummy/Max);

index=find(dummy);

max_relative_error_Var=...
                      max((abs(Var_error(index,1))./Var_feq(index,1)))*100;
max_relative_error_corresponding_absolute_error_Var=...
                                              max(abs(Var_error(index,1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot all relevant data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot 1 has 3 subplots : E0_feq, E0_pmde and E0_error. 
figure(1)
subplot(3,1,1),plot(t,E0_feq);
xlabel('Time');
ylabel('E0(feq)');
subplot(3,1,2),plot(t,E0_pmde);
xlabel('Time');
ylabel('E0(approx)');
subplot(3,1,3),plot(t,E0_error);
xlabel('Time');
ylabel('Error');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input('');
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot 1 has 3 subplots : E0_feq, E0_pmde and E0_error. 
figure(2)
subplot(3,1,1),plot(t,E1_feq);
xlabel('Time');
ylabel('E1(feq)');
subplot(3,1,2),plot(t,E1_pmde);
xlabel('Time');
ylabel('E1(approx)');
subplot(3,1,3),plot(t,E1_error);
xlabel('Time');
ylabel('Error');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input('');
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot 1 has 3 subplots : E0_feq, E0_pmde and E0_error. 
figure(3)
subplot(3,1,1),plot(t,Var_feq);
xlabel('Time');
ylabel('Var(feq)');
subplot(3,1,2),plot(t,Var_pmde);
xlabel('Time');
ylabel('Var (approx)');
subplot(3,1,3),plot(t,Var_error);
xlabel('Time');
ylabel('Error');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following two lines of code were put in by the programmer to clear the
%plots once they have been seen by the user. This has only ornamental
%significance and the user may do away with it. 
input('');
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%