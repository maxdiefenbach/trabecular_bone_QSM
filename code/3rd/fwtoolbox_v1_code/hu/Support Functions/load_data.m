temp1=0;
while temp1 == 0    
   
   [temp_filename,temp_pathname,temp_filterindex]=uigetfile('.mat','SELECT A INPUT DATA FILE');
   fprintf ('\n %s', temp_filename);
   clc;

    switch temp_filename
        %USC
        case 'ObliqueFWPhantom_1ch.mat'
            load ObliqueFWPhantom_1ch.mat; disp ('3T Oblique Phantom Loaded'); temp1=1;
        case 'Ankle_1ch.mat' 
            load Ankle_1ch.mat; disp ('3T Ankle 1ch Loaded');temp1=1;
        case 'Knee_1ch.mat' 
            load Knee_1ch.mat; disp ('3T Knee 1ch Loaded');temp1=1;
        case 'Coronal_Abdomen_3T_15echo.mat'
            load Coronal_Abdomen_3T_15echo.mat; disp ('3T Coronal Abdomen 15echo Loaded');temp1=1;        
        case 'Axial_Liver_3T_15echo.mat'
            load Axial_Liver_3T_15echo.mat; disp ('3T Axial Liver 15echo Loaded');temp1=1;        
        case 'Axial_Liver_3T_15echo(2).mat'
            load Axial_Liver_3T_15echo(2).mat; disp ('3T Axial Liver 15echo Loaded');temp1=1;
        case 'Axial_Abdomen_3T_15echo.mat'
            load Axial_Abdomen_3T_15echo.mat; disp ('3T Axial Abdomen 15echo Loaded');temp1=1;       
        case 'PhantomBottles_3T_3echo.mat'
            load PhantomBottles_3T_3echo.mat; disp ('3T Phantom Bottles Loaded');temp1=1;       
        case 'Bacon_3T_4echo.mat'
            load Bacon_3T_4echo.mat; disp ('3T Bacon 4echo Loaded');temp1=1;   
        
        case 'Ankle_8ch.mat'
            load Ankle_8ch.mat; disp ('3T Ankle 8ch Loaded'); temp1=1;
 
        %Peter Kellman (NIH)
        case 'PKdata1.mat'
            load PKdata1.mat;  disp ('1.5T Short Axis 4echo Loaded');temp1=1;
        case 'PKdata2.mat'
            load PKdata2.mat;  disp ('1.5T 4 Chamber 4echo Loaded'); temp1=1;
        case 'PKdata3.mat'
            load PKdata3.mat;  disp ('1.5T Short Axis 3echo Loaded'); temp1=1;
        case 'PKdata4.mat'
            load PKdata4.mat;  disp ('1.5T 4 Chamber 3echo Loaded'); temp1=1;
        case 'PKdata5.mat'
            load PKdata5.mat;  disp ('1.5T Short Axis 8echo Loaded'); temp1=1;
        case 'PKdata6.mat'
            load PKdata6.mat;  disp ('1.5T 4 Chamber 8echo Loaded'); temp1=1;
        case 'PKdata7.mat'
            load PKdata7.mat;  disp ('1.5T Coronal 3echo Loaded'); temp1=1;
        case 'PKdata8.mat'
            load PKdata8.mat;  disp ('1.5T Coronal 4echo Loaded'); temp1=1;
        case 'PKdata9.mat'
            load PKdata9.mat;  disp ('1.5T Coronal 8echo Loaded'); temp1=1;
        case 'PKdata10.mat'
            load PKdata10.mat; disp ('1.5T Short Axis 3echo DBprep Loaded'); temp1=1;
        case 'PKdata11.mat'
            load PKdata11.mat; disp ('1.5T 4 Chamber 4echo DBprep Loaded'); temp1=1;
        case 'PKdata12.mat'
            load PKdata12.mat; disp ('1.5T Short Axis 3echo DBprep Loaded'); temp1=1;
        case 'PKdata13.mat'
            load PKdata13.mat; disp ('1.5T 4 Chamber 3echo DBprep Loaded'); temp1=1;
        case 'PKdata14.mat'
            load PKdata14.mat; disp ('1.5T 4 Chamber 8echo DBprep Loaded'); temp1=1;
        case 'PKdata15.mat'
            load PKdata15.mat; disp ('1.5T Short Axis 8echo DBprep Loaded'); temp1=1;
            
        %acetone phantom     
        case 'Acetone_1ch.mat'
            load Acetone_1ch.mat; disp ('3T Acetone Phantom'); temp1=1;
            
        otherwise
            disp ('invalid selection');
    end %switch

end %while

clear temp*;
temp = exist ('imDataParams','var');
if (temp == 0)
    clear imDataParams;
    imDataParams = data;
    clear data;
end
clear temp*;