define(["require", "exports"], function (require, exports) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    var DesignUtils;
    (function (DesignUtils) {
        function downloadSVG(svg_element, filename) {
            const s = new XMLSerializer().serializeToString(svg_element);
            const d = "data:image/svg+xml;base64," + window.btoa(s);
            downloadURI(d, filename);
        }
        DesignUtils.downloadSVG = downloadSVG;
        function downloadText(text, filename) {
            const d = "data:text/plain;base64," + window.btoa(text);
            downloadURI(d, filename);
        }
        DesignUtils.downloadText = downloadText;
        function downloadCSV(csv, filename) {
            const d = "data:text/csv;base64," + window.btoa(csv);
            downloadURI(d, filename);
        }
        DesignUtils.downloadCSV = downloadCSV;
        function downloadURI(uri, name) {
            const link = document.createElement("a");
            link.download = name;
            link.href = uri;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
        DesignUtils.downloadURI = downloadURI;
        function downloadPNG(svg_element, filename) {
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
        DesignUtils.downloadPNG = downloadPNG;
    })(DesignUtils = exports.DesignUtils || (exports.DesignUtils = {}));
});
