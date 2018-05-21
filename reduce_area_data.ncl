;================================================================================
;   Copyright 2018 Hannah M. Christensen
;
;   Licensed under the Apache License, Version 2.0 (the "License");
;   you may not use this file except in compliance with the License.
;   You may obtain a copy of the License at
;
;       http://www.apache.org/licenses/LICENSE-2.0
;
;   Unless required by applicable law or agreed to in writing, software
;   distributed under the License is distributed on an "AS IS" BASIS,
;   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;   See the License for the specific language governing permissions and
;   limitations under the License.
;================================================================================


; reduce_area_data.ncl

function reduce_area_data(data_in:float,flag_subset:logical,set_lat_min:float,set_lat_max:float,set_lon_min:float,set_lon_max:float,no_dims:integer)

local    data_in,flag_subset, set_lat_min, set_lat_max, set_lon_min, set_lon_max,no_dims

begin

    if (flag_subset) then
       print("    reduce area data")
       if (no_dims.eq.3) then
          data_out = data_in(time|:,{lat|set_lat_min:set_lat_max},{lon|set_lon_min:set_lon_max})
       else
          if (no_dims.eq.4) then
             data_out = data_in(time|:,eta_hh|:,{lat|set_lat_min:set_lat_max},{lon|set_lon_min:set_lon_max})
          else
             print("WRONG NO DIMS")
          end if
       end if
    else
       data_out = data_in
    end if  

    return(data_out)

end
