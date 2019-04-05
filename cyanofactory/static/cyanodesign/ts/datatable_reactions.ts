import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

$.fn.dataTable();

let table_template = document.createElement('template');
table_template.innerHTML = `
<table id="cyano-enzyme-list" class="table table-striped table-hover">
    <thead>
        <tr>
            <th>Name</th>
            <th>Reaction</th>
            <th>Constraint</th>
            <th>Active</th>
            <th></th>
        </tr>
    </thead>
</table>
`;

let button_template = document.createElement("template");
button_template.innerHTML = `
<div class="btn-group">
    <button type="button" class="btn btn-primary create-enzyme-button">
    Create new reaction
    </button>
    <button type="button" class="btn btn-primary dropdown-toggle" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
        <span class="caret"></span>
        <span class="sr-only">Toggle Dropdown</span>
    </button>
    <ul class="dropdown-menu">
        <li><a href="#" id="create-enzyme-bulk-button">Bulk add reactions</a></li>
    </ul>
</div>
`;

export class Table {
    private tbl: jQuery<HTMLElement>;
    private buttons: jQuery<HTMLElement>;

    constructor(where: HTMLElement) {
        this.tbl = $(document.importNode(table_template,true));
        this.buttons = $(document.importNode(button_template, true));

        $(where).DataTable();
    }
}