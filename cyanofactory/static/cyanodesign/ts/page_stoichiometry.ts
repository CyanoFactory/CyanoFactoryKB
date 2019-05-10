import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

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

export class Page {
    readonly app: app.AppManager;
    readonly datatable: DataTables.Api;
    readonly source_element: HTMLElement;
    readonly table_element: HTMLElement;

    constructor(where: HTMLElement, app: app.AppManager) {
        this.source_element = where;
        where.appendChild(template.content.cloneNode(true));

        this.table_element = <HTMLElement>where.getElementsByClassName("cyano-chemical-list")[0]!;
        this.app = app;
    }

    init() {
        //this.datatable.clear();
        //this.datatable.rows.add(this.model.reactions);
        //this.datatable.draw();
    }
}
