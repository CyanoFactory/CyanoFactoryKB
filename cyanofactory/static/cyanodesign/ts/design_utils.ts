export namespace DesignUtils {
    export function downloadSVG(svg_element: Node, filename: string): void {
        const s = new XMLSerializer().serializeToString(svg_element);
        const d: string = "data:image/svg+xml;base64," + window.btoa(s);
        downloadURI(d, filename);
    }

    export function downloadText(text: string, filename: string): void {
        const d: string = "data:text/plain;base64," + window.btoa(text);
        downloadURI(d, filename);
    }

    export function downloadCSV(csv: string, filename: string): void {
        const d: string = "data:text/csv;base64," + window.btoa(csv);
        downloadURI(d, filename);
    }

    export function downloadJson(json: string, filename: string): void {
        const d: string = "data:application/json;base64," + window.btoa(json);
        downloadURI(d, filename);
    }

    export function downloadURI(uri: string, name: string) {
        const link = document.createElement("a");
        link.download = name;
        link.href = uri;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }

    export function downloadPNG(svg_element: Node, filename: string) {
        const w = svg_element.getBBox().width;
        const h = svg_element.getBBox().height;
        const can = document.createElement('canvas');
        const ctx = can.getContext('2d');
        const loader = new Image;

        loader.width = can.width = w * 4;
        loader.height = can.height = h * 4;
        loader.onload = function () {
            ctx.drawImage(loader, 0, 0, loader.width, loader.height);
            downloadURI(can.toDataURL("image/png"), filename);
        };
        const svgAsXML = (new XMLSerializer).serializeToString(svg_element);
        loader.src = 'data:image/svg+xml,' + encodeURIComponent(svgAsXML);
    }
}

