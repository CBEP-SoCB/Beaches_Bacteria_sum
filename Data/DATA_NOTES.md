# Data on Bacteria LEvels at Casco Bay Beaches

### `beach_locations.csv`
Geospatial data on the location of monitored beaches in Casco Bay. Includes
data for Winslow Park and two "extra" sampling locations at Willard Beach.  These 
samplng locations were seldom sampled over the last few years, so we did not
end up using related data in State of Casco Bay.

Column Name     | Contents                               | Units                         
----------------|----------------------------------------|------
Town            | Name of town or city where beach is located | 
Beach_Name      | Name of beach                          |
SamplePoint     | Sample location codes                  | Alphanumeric
Latitude        | Latitude  WGS 1984                     | Decimal degrees
Longitude        Longitude, WGS 1984                     | Decimal degrees




### `beaches_data.csv`
This file contains the bacteria data from Maine DEP.  The file was significantly
reorganized from teh source data we received from Maine DEP, and the data
provided to us was pulled directly from DEP's "EGAD" data system.  EGAD data is
in "Long" data format with a great deal of metadata provided as additional data
rows, with spcific flags. For details of how we reorganized the data, see the
file `beaches_data_review_and_prep.Rmd` in the related complete data and code
archive.

[here](https://github.com/CBEP-SoCB-Details/Beaches_Bacteria.git).

Column Name     | Contents                               | Units              
----------------|----------------------------------------|-------------------
SiteCode        | Matches "SamplePoint" from locations data |  
sdatetime       | Date and time of sample collection     | yyyy-mm-ddTHH:MM:SS:Z Local clock time (Eastern daylight savings)  
sdate           | Date of sample collection              | mm/dd/yyyy  
Year            | Year of sample collection              | four digit integer  
Month           | Month of sample collection             | integer, 1 to 12  
DOY             | Day of the year                        | 1 - 366  
Sample_ID       | Unique sample identifier, combining location and date |   
Sample_Qualifier| MIssing ("NA"), "Not Applicable" or "Reanalysis" |  
Enterococci     | Bacteria concentration, via "MPN" method | Nominally, Colony Forming Units per 100 ml   
Reporting_Limit | Method reporting limit  (Usually 1 or 10)    | CFU / 100 ml  
Lab_Qualifier   | Various data quality flags.  "U" indicates non-detect. |  
Bacteria        | As Entrococci, but non-detects replaced by reporting limit |   
Censored_Flag   | Flag indicating whether the value in "Bacteria" is the Reporting limit instead of observed value. |   
Rain24          | Rainfall in last 24 hours (Used prior to 2008) | Inches
Rain48          | Rainfall in last 48 hours (Used 2008 to the present) | Inches
Salinity        | Salinity (From sodium)                | Parts per thousand
Air_Temp        | Temperature of the air at time of sample collection  | Degrees Celsius
Water_Temp      | Temperature of the water  at time of sample collection  | Degrees Celsius
Weather         | Narrative of weather at sample collection  | 'CLEAR', 'PARTLY CLOUDY', 'OVERCAST', 'RAIN'  
Past24HR_Weather | Used prior to 2008 | 'HEAVY RAIN', 'MEDIUM RAIN', 'LIGHT RAIN'   
Past48HR_Weather | Used 2008 to the present | 'HEAVY RAIN', 'MEDIUM RAIN', 'LIGHT RAIN', 'NO RAIN'   
Tide_Stage       | 'HIGH', 'HIGH EBB', 'EBB', 'LOW EBB', 'LOW', 'LOW FLOOD', 'FLOOD', 'HIGH FLOOD' |  
Water_Surface   | 'CALM', 'ROUGH'                      |  
Current         | 'SLOW CURRENT', 'MEDIUM CURRENT', 'RAPID CURRENT' |  



