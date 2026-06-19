# -*- coding: utf-8 -*-
"""
Script used to collate results from MC runs for comparison

Created on Wed Jan  8 10:47:25 2025

@author: ChenL
"""

import pandas as pd
import os

# Insert the folder and report list filenames below
folder = r'C:\PythonProjects\fredhaigh\FHD_2024\03_Design\run\E011_12\report'
report_list_files = ['_Report_List.xlsx']
warnings = []


def main():
    for report_list in report_list_files:
        # Import all the data
        filepath = os.path.join(folder, report_list)
        print('Opening report list:', filepath)
        
        report_info = pd.read_excel(filepath, 
                                  sheet_name='reports',
                                  dtype=pd.StringDtype())
        report_info.index = report_info.index.map(str)
        report_info = report_info.loc[report_info['Include'] == 'yes']
        
        urbs_info = pd.read_excel(filepath, 
                                  sheet_name='URBS',
                                  index_col='Key',
                                  dtype=pd.StringDtype())
        urbs_info.index = urbs_info.index.map(str)
        
        # Separate results into spreadsheets
        for report in report_info.iterrows():
            print('\nCreating table', report[0])
            this_report_info = report[1]
            print(this_report_info)
            
            # Get URBS result keys for this report
            urbs_lst = string_to_list(this_report_info['URBS'])
            print('URBS keys:', urbs_lst)
            
            # Create report for results
            create_report(this_report_info, urbs_info.loc[urbs_lst])
        
        ext = '.xlsx'
        if ext in filepath:
            warning_file = filepath.replace('.xlsx', 'errors.txt')
            with open(warning_file, 'w') as f:
                for warning in warnings:
                    f.write(warning+'\n')

def string_to_list(input_str):
    if pd.isna(input_str):
        return None
    
    if ',' in str(input_str):
        input_lst = [n.strip() for n in input_str.split(',')]
    else:
        input_lst = [input_str]
    return input_lst

def collate_urbs_results(filepathtemplate, result_types, do_output=True):
    print('URBS result types:', result_types)
    outpath = filepathtemplate.replace('~type~.csv', 'report.xlsx')
    writer = pd.ExcelWriter(outpath, engine='openpyxl', mode='w')
    for result_type in result_types:
        filepath = filepathtemplate.replace('~type~', result_type)
        print('Opening URBS result file:', filepath)
        try:
            df = pd.read_csv(filepath)
            df.set_index('aep (1 in x)', inplace=True)
            if result_type == 'inflow':
                df_inflow = df.copy()
            elif result_type == 'level':
                df_level = df[['max','critical_duration']]
            else:
                df_outflow = df[['max']]
            df.to_excel(writer, sheet_name=result_type)
        except Exception:
            warning = f'WARNING: failed to open: {filepath}'
            warnings.append(warning)
            print(warning)
    
    print('\nCollating results for report')
    col_names = ['Peak inflow (m3/s)','Inflow critical duration (h)','Peak inflow* (m3/s)', 'Peak outfow (m3/s)', 
                 'Peak level (m AHD)', 'Level critical duration (h)']
    df_report = pd.DataFrame(index=df_inflow.index, columns=col_names)
    df_report[['Peak inflow (m3/s)']] = df_inflow[['max']]
    df_report['Inflow critical duration (h)'] = df_inflow.critical_duration
    df_report[['Peak outfow (m3/s)']] = df_outflow[['max']]
    df_report[['Peak level (m AHD)']] = df_level[['max']]
    df_report['Level critical duration (h)'] = df_level.critical_duration
    df_report['Peak inflow* (m3/s)'] = [df_inflow[df_report['Level critical duration (h)'][i]][i] for i in df.index] # Find peak inflow for level critical duration
    df_report.to_excel(writer, sheet_name='report')
    
    writer.close()
    
    return df_report

def create_report(report_info, urbs_info):
    outpath =  os.path.join(report_info['Folder'], report_info['Filename'])
    outpath = f'{outpath}.xlsx'
    writer = pd.ExcelWriter(outpath, engine='openpyxl', mode='w')
    
    level = pd.DataFrame()
    outflow = pd.DataFrame()
    critical_dur = pd.DataFrame()
    inflow = pd.DataFrame()
    inflow_cd = pd.DataFrame()
    
    for this_urbs_info in urbs_info.iterrows():
        this_urbs = this_urbs_info[1]
        print('\nURBS simulation info:', this_urbs_info[0])
        print(this_urbs)
        filepath = os.path.join(this_urbs['Folder'], this_urbs['Filename'])
        types = string_to_list(this_urbs['Type'])
        urbs_results = collate_urbs_results(filepath,types, do_output=True)
        
        if urbs_results is not None:
            urbs_results.to_excel(writer,sheet_name=this_urbs['Label'])
            
            # Combine results
            level = pd.concat([level,urbs_results['Peak level (m AHD)']],axis=1)
            level.columns = [*level.columns[:-1], this_urbs['Label']]
            
            outflow = pd.concat([outflow,urbs_results['Peak outfow (m3/s)']],axis=1)
            outflow.columns = [*outflow.columns[:-1], this_urbs['Label']]

            critical_dur = pd.concat([critical_dur,urbs_results['Level critical duration (h)']],axis=1)
            critical_dur.columns = [*critical_dur.columns[:-1], this_urbs['Label']]
            
            inflow = pd.concat([inflow,urbs_results['Peak inflow (m3/s)']],axis=1)
            inflow.columns = [*inflow.columns[:-1], this_urbs['Label']]
            
            inflow_cd = pd.concat([inflow_cd,urbs_results['Inflow critical duration (h)']],axis=1)
            inflow_cd.columns = [*inflow_cd.columns[:-1], this_urbs['Label']]

    level.index.name = 'Level'            
    level.to_excel(writer,sheet_name='comparison', startrow=0)
    
    outflow.index.name = 'Outflow'
    outflow.to_excel(writer,sheet_name='comparison', startrow=len(outflow)+2)
    
    critical_dur.index.name = 'Level CD'
    critical_dur.to_excel(writer,sheet_name='comparison', startrow=(len(outflow)+2)*2)
    
    inflow.index.name = 'Inflow'
    inflow.to_excel(writer,sheet_name='comparison', startrow=(len(outflow)+2)*3)

    inflow_cd.index.name = 'Inflow CD'
    inflow_cd.to_excel(writer,sheet_name='comparison', startrow=(len(outflow)+2)*4)
    
    writer.close()
    
    print('Excel closed')

if __name__ == "__main__":
    main()