import plotly.express as px

def drawPicture(dataframe,col_name, row_name,colorattribute,save_file,celltype_colors =  ("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A80",
    "#377F1C", "#4DAFAE", "#984F07", "#FF7F64", "#FFFF97", "#A6568C", "#F78223", "#9999FD", "#E41AE4", "#377F80", "#4DB012",
    "#984F6B", "#FF7FC8",  "#A656F0", "#F78287", "#999A61", "#E41B48", "#377FE4", "#4DB076", "#984FCF", "#FF802C",
    "#00005F", "#A65754", "#F782EB", "#999AC5", "#E41BAC", "#378048", "#4DB0DA", "#985033", "#FF8090", "#0000C3", "#A657B8",
    "#F7834F", "#999B29", "#E41C10", "#3780AC", "#4DB13E", "#985097", "#FF80F4", "#000127", "#A6581C", "#F783B3", "#999B8D",
    "#E41C74", "#378110", "#4DB1A2", "#9850FB", "#FF8158", "#00018B", "#A65880", "#F78417", "#999BF1", "#E41CD8", "#378174",
    "#4DB206", "#98515F", "#FF81BC", "#0001EF", "#A658E4", "#F7847B", "#999C55", "#E41D3C", "#3781D8", "#4DB26A", "#9851C3",
    "#FF8220", "#000253", "#A65948", "#F784DF", "#999CB9", "#E41DA0", "#37823C", "#4DB2CE", "#985227", "#FF8284", "#0002B7",
    "#A659AC", "#F78543", "#999D1D", "#E41E04", "#3782A0", "#4DB332", "#98528B", "#FF82E8", "#00031B", "#A65A10", "#F785A7",
    "#999D81", "#E41E68", "#378304", "#4DB396", "#9852EF", "#FF834C", "#00037F", "#A65A74", "#F7860B", "#999DE5", "#E41ECC",
    "#378368", "#4DB3FA", "#985353", "#FF83B0", "#0003E3", "#A65AD8", "#F7866F", "#999E49", "#E41F30", "#3783CC", "#4DB45E",
    "#9853B7", "#FF8414", "#000447", "#A65B3C", "#F786D3", "#999EAD", "#E41F94", "#378430", "#4DB4C2", "#98541B", "#FF8478",
    "#0004AB", "#A65BA0", "#F78737", "#999F11", "#E41FF8", "#378494", "#4DB526", "#98547F", "#FF84DC", "#00050F", "#A65C04",
    "#F7879B", "#999F75"
    ),width = 1000,height = 1000,is_show = True,is_save = False,save_type = "pdf",point_size = 3):
    size = len(set(dataframe[colorattribute]))
    dataframe.sort_values(by = colorattribute,inplace=True, ascending=True)
    length_col = max(dataframe[col_name]) - min(dataframe[col_name])
    length_row = max(dataframe[row_name]) - min(dataframe[row_name])
    max_length = max(length_col,length_row) + 2
    fig = px.scatter(dataframe, x = col_name, y= row_name,color = colorattribute,color_discrete_sequence=celltype_colors)
    fig.update_traces(marker_size = point_size)
    fig.update_layout(
        xaxis = dict(
            tickmode = 'linear',
            tick0 = min(dataframe[row_name]), 
            dtick = max_length 
        ),
        yaxis = dict(
            tickmode = 'linear',
            tick0 = min(dataframe[col_name]),
            dtick = max_length
        ),

    )
    fig.update_layout(
        autosize=False,
        width=width,
        height = height)
    if(is_show):
        fig.show()
    if(is_save):
        if(save_type == "png"):
            fig.write_image(save_file)
        if(save_type == "pdf"):
            fig.write_image(save_file)
        if(save_type == "html"):
            fig.write_html(save_file)
