import dash_html_components as html
import dash_core_components as dcc

def Header():
    return html.Div([
        get_logo()
    ])

def get_logo():
    logo = html.Div([

        html.Div([
            html.Img(src='/assets/DetCat.png', id='logo', height='120', width='100')
        ], style={
                    'width': '10%',
                    'display': 'inline-block',
                    'align':'middle'

                }, className="w3-bar-item w3-button w3-blue"),
        html.Div([
            html.H1('DetCat: Detecting Categorical Outliers'),
            html.H1(' in Relational Datasets'),
            
        ], style={
                    'white-space':'pre',
                    'height': '120',
                    'padding-top':'1em',
                    'display': 'inline-block',
                    'vertical-align':'center'
                })

    ], className="w3-bar w3-blue")


    # ], className="row gs-header")
    return logo

def Footer():
    # return html.Div([
    #         html.Img(src='/assets/Logo_red.png', id='heart', height='120', width='165')
    #     ], style={
    #                 'width': '15%',
    #                 'display': 'inline-block',
    #                 'align':'middle'

    #             }, className="w3-bar-item w3-button w3-white"),
        html.Div([
            html.A('Data Intensive Systems - UU', href='https://www.uu.nl/en/research/ai-data-science/data-intensive-systems')
        ])
