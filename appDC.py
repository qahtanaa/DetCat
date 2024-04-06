import base64
import datetime
import io
import os 
from os import listdir
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import ntpath
import csv
import pandas as pd
# from extra import find_pfds_csv
from detective_cat import run_dc
# from dash import Dash, Input, Output, callback

from components import Header

DATA_FOLDER = "./Datasets/Dirty_Data"


# important link https://dash.plot.ly/datatable/interactivity
# external CSS stylesheets
external_stylesheets = [
    'https://codepen.io/chriddyp/pen/bWLwgP.css',
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO',
        'crossorigin': 'anonymous'
    }
]



app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.config['suppress_callback_exceptions']=True
app.title = 'DetCat'

gdf = None
gresults = None
param_dict = dict()
param_dict["results_main_dir"] = "../Results/"
sim_functions = ['Levenshtein distance', 'Cust. Lev. distance']



def read_table(tab_name):
    t_name = ntpath.basename(tab_name)
    try:
        df = pd.read_csv(filepath_or_buffer=tab_name, dtype=object, delimiter=',', low_memory=False,
                         quoting=csv.QUOTE_ALL, doublequote=True)
    except ValueError:
        try:
            df = pd.read_csv(filepath_or_buffer=tab_name, dtype=object, delimiter=',', low_memory=False,
                             quoting=csv.QUOTE_ALL, doublequote=True, encoding="ISO-8859-1")
        except:
            print("Error reading csv file .. file encoding is not recognizable")
            return None
    return df

def get_csv_files(folder):
    csv_files = []
    folder_contents = listdir(folder)
    for item in folder_contents:
        if item.endswith('.csv'):
            csv_files.append(item)
    return csv_files


def dynamic_page():
    global gdf
    column_options = list('All')

    return html.Div([
        Header(),
        html.Div([    
            html.Div([
                dcc.Dropdown(
                        id='uploaded-datasets',
                        options=[{'label': i, 'value': i} for i in get_csv_files(DATA_FOLDER)],
                        placeholder='Select a Dataset',
                    )
                ],style={
                        'width': '300px',
                        'display': 'inline-block',
                        'margin-left': '25px',
                        'font-size': '150%'
                    } 
                
            ), 
            html.Div([
                dcc.Upload(
                    id='upload-data',
                    children = html.Div([
                        html.Button(' Upload ', className='fa fa-upload')
                    ],style={
                        'backgroundColor':'green',
                        'color':'white',
                        'margin-left': '5px',
                        'font-size': '150%'
                    }),
                    
                    multiple=False
                ),
            ]),
        ], className="row",
            style={
                'width': '100%',
                'height':'50px',
                'borderWidth': '1px',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin-left': '25px',
                'margin-top': '10px',
            }),
        html.Hr(),  # horizontal line
        html.Div(id='output-data-upload'),
        html.Div(id='output-data-dropdown',
            style={
                'width': '100%',
                'height': '440px',
                'borderWidth': '1px',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin-left': '50px',
                'margin-right': '25px',
                'overflowY': 'scroll',
            }),
        html.Hr(),  # horizontal line        
        html.Div([
            html.H2('Parameter Settings',
                    style={'display':'inline-block','margin-left':20, 
                        'float':"left", 'font-size': '250%',
                'font-family': 'Times New Roman'}),
            html.Div([
                html.B('Min Coverage',
                    style={'display':'inline-block','margin-left':80, 
                        'color':'green', 'float':"left"}),
                dcc.Input(id='MCoverage',
                        type='text',
                        placeholder='99', 
                        value='99',
                        style={'display':'inline-block', 'margin-left':105, 'float':"left",
                        'border': '1px solid black'}),
                ], style={'display':'inline-block', 
                        'margin-left': '5px', 
                        'width': '100%', 
                        'float':"left", 
                        'font-size': '150%'}),                       

            html.Div([
                html.B('Max Length',
                    style={'display':'inline-block','margin-left':80, 
                        'color':'green', 'float':"left"}),
                dcc.Input(id='MLen',
                        type='text',
                        placeholder='10', 
                        value='10',
                        style={'display':'inline-block', 'margin-left':130, 'float':"left",
                        'border': '1px solid black'}),
                ], style={'display':'inline-block', 
                        'margin-left': '5px', 
                        'width': '100%', 
                        'float':"left", 
                        'font-size': '150%'}),
            
            html.Div([
                html.B('Min Cov. per Pattern',
                    style={'display':'inline-block','margin-left':80, 
                        'color':'green', 'float':"left"}),
                dcc.Input(id='MCDP',
                        type='text',
                        placeholder='0.1', 
                        value='0.1',
                        style={'display':'inline-block', 'margin-left':20, 'float':"left",
                        'border': '1px solid black'}),
                ], style={'display':'inline-block', 
                        'margin-left': '5px', 
                        'width': '100%', 
                        'float':"left", 
                        'font-size': '150%'}),
            
            html.Div([
                html.B('Distance Function',
                    style={'display':'inline-block','margin-left':80, 
                        'color':'green', 'float':"left", 'font-size': '150%'}),
                dcc.Dropdown(
                        id='sim_function',
                        options=[{'label': i, 'value': i} for i in sim_functions],
                        # placeholder='Distance Function',
                        style={'display':'inline-block', 'margin-left':25, 'float':"left",
                        'width': '300px', 'height':'40px', 'font-size': '130%', 'vertical-align':'middle'}
                    ),
                ],style={'display':'inline-block', 
                        'margin-left': '5px', 
                        'width': '100%', 
                        'float':"left"}),
            
            html.Div([
                html.Button('Run DetCat', className='fa', id='button', 
                    style={
                            'backgroundColor':'green',
                            'color':'white',
                            'width': '200',
                            'flow':'right',
                            'margin-left': '200px',
                            'float':'left',
                            'font-size': '200%'                            
                        })
                ], style={'width':'100%'}),
            html.Hr(),
        ], className="row",
            style={
                'width': '100%',
                'height':'300px',
                'borderWidth': '1px',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin-left': '25px',
                'margin-top': '10px',
            }),
        
        html.Hr(),  # horizontal line
        html.Div(id = 'output-results',
            style={
                'width': '100%',
                'height': '200px',
                'borderWidth': '1px',
                'borderRadius': '5px',
                'textAlign': 'left',
                'margin-left': '55px',
                'margin-right': '25px',
                'float':'right',
                'margin-top': '40px'
            }),
        # Footer()
    ], className='body')

app.layout = dynamic_page



def upload_contents(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    new_file_name = ''
    try:
        if filename.endswith('csv'):
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
            if filename.endswith('xls'):
                new_file_name = filename.replace('.xls', '.csv')
            if filename.endswith('xlsx'):
                new_file_name = filename.replace('.xlsx', '.csv')
        filename2 = os.path.join(DATA_FOLDER, filename)
        df.to_csv(filename2, sep=',', encoding='latin', index=False, quoting=csv.QUOTE_ALL, doublequote=True)
    except Exception as e:
        print(e)
        return html.Div([
            ''
        ])
    return html.Div([
            ''
        ])


def parse_contents(filename):
    global gdf
    try:
        if filename.endswith('csv'):
            # Assume that the user uploaded a CSV file
            filename2 = os.path.join(DATA_FOLDER, filename)
            df = read_table(filename2)
            gdf = df
        # elif 'xls' in filename:
        #     # Assume that the user uploaded an excel file
        #     df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns],
            css=[{
                'selector': '.dash-cell div.dash-cell-value',
                'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
            }],
            style_cell={
                'whiteSpace': 'no-wrap',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'maxWidth': 0,
                'font_size': '200%',
                'textAlign':'left'
            },
            style_cell_conditional=[{
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }],
            style_header={
                'backgroundColor': 'white',
                'fontWeight': 'bold'
            },
            style_table={
                'max_rows_in_viewport':15,
                'maxHeight': '400px',
                'overflowY': 'scroll'
            },
            pagination_settings={
            "current_page": 0,
            "page_size": 50,
            },
        ),
        # html.Hr(),  # horizontal line
    ], className='ui-grid-resize-columns ui-grid-pagination')


@app.callback([Output('output-data-upload', 'children'),
                Output('uploaded-datasets', 'options'),
                Output('uploaded-datasets', 'value')],
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def update_output_data(content, fname, modified):
    if content:
        grid = upload_contents(content, fname)
        print(content[0])
        options=[{'label': i, 'value': i} for i in get_csv_files(DATA_FOLDER)]
        value=fname
        return grid, options, value
    else:
        options=[{'label': i, 'value': i} for i in get_csv_files(DATA_FOLDER)]
        return html.Div(['']), options, ''


@app.callback(Output('output-data-dropdown', 'children'),
              [Input('uploaded-datasets', 'value')])
def output_dropdown(fname):
    options=[{'label': i, 'value': i} for i in get_csv_files(DATA_FOLDER)]
    if fname:
        grid = parse_contents(fname)
        return grid
    else:
        return html.Div([''])
        
@app.callback(
    Output('output-results', 'children'),
    [Input('button', 'n_clicks')],
    [State('uploaded-datasets', 'value'),
    State('MCoverage', 'value'),
    State('MLen', 'value'),
    State('MCDP', 'value'),
    State('sim_function', 'value')])
def update_output_discovery(n_clicks, fname, min_coverage, MLen, MCDP, sim_function):
    
    global gresults

    if fname:
        param_dict["tab_name"] = os.path.join(DATA_FOLDER, fname)
        pmcov = float(min_coverage) / 100.0 if 0 <= float(min_coverage) <= 100 else 0.95
        param_dict["min_acceptable_coverage"] = pmcov 
        pmlen = float(MLen) if 5 <= float(MLen) <= 20 else 15
        param_dict["MaxLen"] = pmlen
        pMCDP = float(MCDP) if 0.01 <= float(MCDP) <= 0.5 else 0.1
        param_dict["MCDP"] = pMCDP
        param_dict["sim_function"] = sim_function
        gresults = run_dc(param_dict)
        if len(gresults) > 0:
            columns = gresults.keys()
            data = list()
            for kk in range(len(data)):
                data.remove(data[0])  # make sure the data list is empty  
            for c in columns:
                data.append((c, gresults[c]["dType"])) 
            if len(columns) > 0:
                init_df = pd.DataFrame(data, columns = ['columns', 'Type'])
                return html.Div([
                    dash_table.DataTable(
                    data = init_df.to_dict('records'),
                    columns=[{'name': i, 'id': i} for i in init_df.columns],
                    id='syn-table',
                    row_selectable="single",
                    # row_deletable=True,
                    css=[{
                        'selector': '.dash-cell div.dash-cell-value',
                        'rule': 'display: inline; white-space: inherit; margin-left: 20px; overflow: inherit; text-overflow: inherit;'
                    }],
                    style_cell={
                        'whiteSpace': 'no-wrap',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'maxWidth': '500px',
                        'margin-left':'55px',
                        'textAlign':'left',
                        'font-size': '200%',
                    },
                    style_cell_conditional=[{
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(248, 248, 248)'
                    }],
                    style_header={
                        'backgroundColor': 'white',
                        'fontWeight': 'bold'
                    },
                    style_table={
                        'max_rows_in_viewport':15,
                        'maxHeight': '400px',
                        'maxWidth':'800px',
                        'overflowY': 'scroll',
                        'margin-left': '20px',
                        # 'border': 'thin lightgrey solid',
                    }
                ),
                    # html.H3('Select a dependency to see its tableau of PFDs'),
                    html.Div(id='syn-container', className="six columns"),
                    # html.Div(id='pfds-container-hidden', style={'display':'none'}),
                ], className="row ")      
        else:
            html.Div([''])
    else:
        html.Div(['The data file is missing'])
    if n_clicks:
        return html.Div(['Something goes wrong after {}'.format(n_clicks)])
    return html.Div([''])





@app.callback(
    Output('syn-container', "children"),
    [Input('syn-table', "derived_virtual_data"),
     Input('syn-table', "derived_virtual_selected_rows")])
def update_graphs_patterns(rows, derived_virtual_selected_rows):
    # When the table is first rendered, `derived_virtual_data` and
    # `derived_virtual_selected_rows` will be `None`. This is due to an
    # idiosyncracy in Dash (unsupplied properties are always None and Dash
    # calls the dependent callbacks when the component is first rendered).
    # So, if `rows` is `None`, then the component was just rendered
    # and its value will be the same as the component's dataframe.
    # Instead of setting `None` in here, you could also set
    # `derived_virtual_data=df.to_rows('dict')` when you initialize
    # the component.
    global gresults, gdf
    if not(derived_virtual_selected_rows):
        derived_virtual_selected_rows = []
        return html.Div([
            html.H4('') 
        ], className="six columns")
    else:
        cur_att = rows[derived_virtual_selected_rows[0]]['columns']
        ptrns = gresults[cur_att]['dominant']
        att_ptrns = []
        for kk in range(len(att_ptrns)):
            att_ptrns.remove(att_ptrns[0])  # make sure the list is empty 
        for k in ptrns:
            att_ptrns.append((k, ptrns[k]))
        # print (att_ptrns)
        
        outliers = gresults[cur_att]['outlier patterns']
        att_outliers = []
        for kk in range(len(att_outliers)):
            att_outliers.remove(att_outliers[0])  # make sure the list is emplty 
        for k in outliers:
            att_outliers.append((k, outliers[k]))
        # print (att_ptrns)

        patt_df = pd.DataFrame(att_ptrns, columns=['pattern (p)', '# values represented by p'])
        outliers_df = pd.DataFrame(att_outliers, 
                        columns=['Outlier', 'Frequency'])
        return html.Div([
            dash_table.DataTable(
            data=patt_df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in patt_df.columns],
            id='ptrns-syn-table',
            css=[{
                'selector': '.dash-cell div.dash-cell-value',
                'rule': 'display: inline; white-space: inherit; margin-left: 20px; overflow: inherit; text-overflow: inherit;'
            }],
            style_cell={
                'whiteSpace': 'no-wrap',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'maxWidth': '600px',
                'textAlign':'left',
                'font-size': '200%',
                'font-family': 'Times New Roman'
            },
            style_cell_conditional=[{
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }],
            style_header={
                'backgroundColor': 'white',
                'fontWeight': 'bold'
            },
            style_table={
                'max_rows_in_viewport':15,
                'maxHeight': '400px',
                'maxWidth':'800px',
                'overflowY': 'scroll',
                'margin-left': '20px',
                # 'border': 'thin lightgrey solid',
            }
        ),
            dash_table.DataTable(
            data=outliers_df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in outliers_df.columns],
            id='outliers-syn-table',
            css=[{
                'selector': '.dash-cell div.dash-cell-value',
                'rule': 'display: inline; white-space: inherit; margin-left: 20px; overflow: inherit; text-overflow: inherit;'
            }],
            style_cell={
                'whiteSpace': 'no-wrap',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'maxWidth': '600px',
                'textAlign':'left',
                'font-size': '200%',
                'font-family': 'Times New Roman'
            },
            style_cell_conditional=[{
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }],
            style_header={
                'backgroundColor': 'white',
                'fontWeight': 'bold'
            },
            style_table={
                'max_rows_in_viewport':15,
                'maxHeight': '400px',
                'maxWidth':'800px',
                'overflowY': 'scroll',
                'margin-left': '20px',
                # 'border': 'thin lightgrey solid',
            }
        ),
        ])
            # else:
            #     text = '(' + gdf.columns[derived_virtual_selected_rows[0]] + ') has been ignored because it represents ' 
            #     text += 'a numerical quantity '
            #     return html.Div([
            #         # html.H3('The selected attribute is ( ' + derived_virtual_selected_rows[0] + ' )')
            #         html.H3(text)
            #     ], className="six columns")
    return html.Div([
        # html.H3('The selected attribute is ( ' + derived_virtual_selected_rows[0] + ' )')
        html.H3('')
    ], className="six columns")






# # # # # # # # #
# detail the way that external_css and external_js work and link to alternative method locally hosted
# # # # # # # # #
external_css = ["https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
                "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                "//fonts.googleapis.com/css?family=Raleway:400,300,600",
                "https://codepen.io/bcd/pen/KQrXdb.css",
                "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
                "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css",
                "https://www.w3schools.com/w3css/4/w3.css",
                "/assets/style.css"]

for css in external_css:
    app.css.append_css({"external_url": css})

external_js = ["https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js",
               "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js",
               ]

for js in external_js:
    app.scripts.append_script({"external_url": js})






if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')


