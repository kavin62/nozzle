# bell-nozzle-contour

This project aim is to generate bell nozzle contour in your respective designing software's (like Catia, Ansys, etc).

BEFORE USING THIS CODE PLEASE DOWNLOAD THE 2 EXCEL FILES

1) 'NOZZLEPARA' CONSISTS OF FEW NOZZLE PARAMETERS, MAKE RESPECTIVE CHANGES IN IT TO GET YOUR DESIRED NOZZLE CONTOUR.

    --> AFTER DOWNLOADING THE EXCEL SHEET IN YOUR SYSTEM, COPY THE FILE PATH AND REPLACE IT WITH THE PATH ALREADY PRESENT INSIDE THE CODE(INSIDE READ EXCEL FUNCTION).
    --> df = read_excel(r' ENTER YOUR FILE PATH HERE!! ')
    
2) 'NOZZLEPTS' FILE CONTAINS ALL THE COORDINATES GENERATED FROM THE CODE.
    
    --> BASICALLY AFTER YOU RUN THE CODE, THE COORDINATES GENERATED FROM THE CODE WILL BE STORED IN YOUR DESIRED EXCEL SHEET.
    --> TO STORE IN YOUR DESIRED EXCEL SHEET JUST REPLACE THE FILE NAME 'nozzlepts.xls' with your desired one near the last line of code (ex: wb.save('example.xls') ).
    
NOTE:
Before you run this code , you have to make sure that following python packages are installed:
1.Pandas
2.scipy
3.xlwt
4.matplotlib

