classdef script_manager < handle
    %SCRIPT_MANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties
        output_file_name
    end
    
    properties(SetAccess = private, GetAccess = private)
        output_file_ID
        
    end
    
    methods
        function obj = script_manager(title,description,version_number,date_created)
            % Setup enviroment.
            format shortEng
            % Setup output file.
            obj.output_file_name = strcat(regexprep(title,' ','_'),'.htm');
            obj.output_file_ID = fopen(obj.output_file_name,'w','n','Shift_JIS');
            if obj.output_file_ID == -1;error('Output File Failed to Open.');end;
            obj.writeline_raw('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
            obj.writeline_raw('<html xmlns="http://www.w3.org/1999/xhtml">');
            obj.writeline_raw('<head>');
            obj.writeline_raw(strcat('<title>',title,'</title>'));
            obj.writeline_raw('<meta http-equiv="Content-Language" content="en-us" />');
            obj.writeline_raw('<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />');
            % Use inline CSS.
            obj.writeline_raw('<STYLE>');
            obj.writeline_raw('<!--');
            % Elements.
            obj.writeline_raw('h1 { font-family: Verdana; font-size: 24px; color: #0066FF }');
            obj.writeline_raw('h2 { font-family: Verdana; font-size: 18px; color: #333333 }');
            obj.writeline_raw('body { font-family: Calibri; font-size: 14px; color: #000000 }');
            % Classes.
            obj.writeline_raw('.table_style_1 { width: 100%% }');
            obj.writeline_raw('.table_header_text { font-family: Calibri; font-size: 14px; font-weight: bold; color: #000000 }');
            
            obj.writeline_raw('-->');
            obj.writeline_raw('</STYLE>');
            obj.writeline_raw('</head>');
            obj.writeline_raw('<body>');
            
            % Print Title Block to output file.
            obj.h1(title);
            obj.h2(description);
            obj.writeline('By Dallas Adams - S2837800');
            obj.writeline(strcat('Version : ',version_number));
            obj.writeline(strcat('Created : ',date_created));
            
        end
        
        function make_table(obj,data,column_headers)
            % Open Table.
            obj.writeline_raw('<table class="table_style_1">');
            if nargin > 1
                % Create Header.
                [no_of_headers,null] = size(column_headers);
                header_line = '<tr>';
                for i = 1:no_of_headers
                    header_line = strcat(header_line,'<td class="table_header_text">',char(column_headers(i)),'</td>');
                end
                obj.writeline_raw( strcat(header_line,'</tr>'));
            end
            % Create Table Rows.
            [no_of_rows,no_of_columns] = size(data);
            for j = 1:no_of_rows
                row_line = '<tr>';
                for i = 1:no_of_columns
                    row_line = strcat(row_line,'<td>',eng(data(j,i)),'</td>');
                end
                obj.writeline_raw(strcat(row_line,'</tr>'));
            end
            obj.writeline_raw('</table>');
        end
        
        function h1(obj,title_str)
            fprintf(obj.output_file_ID,'<h1>');
            fprintf(obj.output_file_ID,title_str);
            fprintf(obj.output_file_ID,'</h1>\n');
        end
        
        function h2(obj,title_str)
            fprintf(obj.output_file_ID,'<h2>');
            fprintf(obj.output_file_ID,title_str);
            fprintf(obj.output_file_ID,'</h2>\n');
        end
        
        function writeline(obj,title_str)
            fprintf(obj.output_file_ID,title_str);
            fprintf(obj.output_file_ID,'</br>\n');
        end
        
        function delete(obj)
            % This is the deconstructor. Closes the object.
            obj.writeline_raw('</body>');
            obj.writeline_raw('</html> ');
            fclose(obj.output_file_ID);
            winopen(obj.output_file_name);
            % [status,cmdout] = system(obj.output_file_name);
        end
        
        function writeline_raw(obj,new_line)
            fprintf(obj.output_file_ID,new_line);
            fprintf(obj.output_file_ID,'\n');
        end
        
        
        function error_test(obj)
            try
                
                throw('error 1')
                
            catch ME
                fprintf('ERROR\n');
                fprintf('Identifier:  %s\n',ME.identifier);
                fprintf('Message   :  %s\n',ME.message);
                fprintf('Full Path : %s\n',mfilename('fullpath'));
                fprintf('Class     : %s\n',mfilename('class'));
                
                
                fprintf('Variables :\n');
                whos;
                
                %fprintf('dbstack   : %s\n',dbstack);
                
                %Open the workspace to view.
                workspace
                %throw the error again.
                rethrow(ME);
            end
        end
        
    end
end

