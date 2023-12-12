"""PDF report generation utility methods."""
from pathlib import Path

import pdfkit


def get_table_html(html_table: str, table_name: str) -> str:
    """Return HTML-rendered content with the provided HTML table."""
    return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                h2 {{text-align: center; padding: 10px;}}
                table {{margin: 0 auto; border: 1px solid black; border-collapse: collapse; text-align: center;}}
                th {{font-size: 12pt; padding: 5px; background: #cccccc;}}
                td {{font-size: 10pt; padding: 5px;}}
                tr {{page-break-inside: avoid;}}
                tr:nth-child(even) {{background: #eeeeee;}}
            </style>
        </head>
        <body>
            <h2>{table_name}</h2>
            {html_table}
        </body>
        </html>
    """


def get_image_html(image_path: Path, image_name: str) -> str:
    """Return HTML-rendered content with the provided image."""
    return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                h1 {{text-align: center;}}
                img {{max-height: 800px;}}
                div {{
                    height: 800px;
                    display: -webkit-box;
                    -webkit-box-pack: center; /* Center horizontally */
                    -webkit-box-align: center; /* Center vertically */
                    position: relative;
                }}
            </style>
        </head>
        <body>
            <h1>{image_name}</h1>
            <div>
                <img src="{image_path.as_posix()}">
            </div>
        </body>
        </html>
    """


def html_to_pdf(
    html_string: str,
    pdf_path: str,
    orientation: str = "landscape",
    margin_top: str = "1.5cm",
    margin_bottom: str = "1cm",
    margin_left: str = "1cm",
    margin_right: str = "1cm",
    zoom: int = 1,
) -> None:
    """Create a PDF file from the content of an HTML string."""
    pdfkit.from_string(
        input=html_string,
        output_path=pdf_path,
        options={
            "page-size": "A4",
            "encoding": "UTF-8",
            "orientation": orientation,
            "zoom": zoom,
            "margin-top": margin_top,
            "margin-bottom": margin_bottom,
            "margin-left": margin_left,
            "margin-right": margin_right,
            "enable-local-file-access": None,
        },
    )
