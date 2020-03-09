# import IPython.core.display
import matplotlib as mpl

def clean():
    # set mpl defaults for nice display
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['figure.figsize'] = (18, 6)
    mpl.rcParams['lines.linewidth'] = 1

    # return IPython.core.display.HTML("""
    # <style type="text/css">
    # div.input {
    # width: 105ex; /* about 80 chars + buffer */
    # }
    # div.text_cell {
    # width: 105ex /* instead of 100%, */
    # }
    # div.text_cell_render {
    # font-family: Cambria, Candara, serif;
    # font-size: 14pt
    # line-height: 145%; /* added for some line spacing of text. */
    # width: 105ex; /* instead of 'inherit' for shorter lines */
    # }
    # /* Set the size of the headers */
    # div.text_cell_render h1 {
    # font-size: 20pt;
    # }
    # div.text_cell_render h2 {
    # font-size: 18pt;
    # }
    # .CodeMirror {
    # font-family: Consolas, monospace;
    # width: 105ex;
    # }
    # .rendered_html ol {list-style:decimal; margin: 1em 2em;}
    # </style>""")
