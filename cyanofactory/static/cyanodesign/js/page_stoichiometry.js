define(["require", "exports", "datatables.net"], function (require, exports) {
    "use strict";
    Object.defineProperty(exports, "__esModule", { value: true });
    let template = document.createElement('template');
    template.innerHTML = `
<table class="cyano-chemical-list table table-striped table-hover">
    <thead>
        <tr>
            <th>Left</th>
            <th>Right</th>
            <th>Difference</th>
        </tr>
    </thead>
</table>
`;
    class Page {
        constructor(where, app) {
            this.source_element = where;
            where.appendChild(template.content.cloneNode(true));
            this.table_element = where.getElementsByClassName("cyano-chemical-list")[0];
            this.app = app;
        }
        update() {
            //this.datatable.clear();
            //this.datatable.rows.add(this.model.reactions);
            //this.datatable.draw();
        }
    }
    exports.Page = Page;
});
