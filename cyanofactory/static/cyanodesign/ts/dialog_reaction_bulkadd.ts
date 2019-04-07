import * as mm from "./metabolic_model";
import * as $ from "jquery";
import "datatables.net";

let template = document.createElement('template');
template.innerHTML = `
<div class="dialog-reaction-bulkadd modal fade" tabindex="-1" role="dialog" aria-labelledby="dialog-reaction-bulkadd-label"
     aria-hidden="true">
    <div class="modal-dialog">
        <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span
                        aria-hidden="true">&times;</span></button>
                <h4 class="modal-title" id="dialog-reaction-bulkadd-label">Bulk add reactions</h4>
            </div>
            <div class="modal-body">
                <form>
                    <fieldset>
                        Enter reactions in BioOpt format for adding them to the model. Constraints can be put after the reactions.<br>
                        <div class="form-group">
                            <label class="control-label" for="enzyme-bulkdata">BioOpt data:</label>
                            <textarea class="form-control" id="enzyme-bulkdata" placeholder="reac1: a + 2 b -> 3.5 c [-1, 1]" rows="5"></textarea>
                        </div>
                        <div class="form-group" id="enzyme-bulkadd-preview">

                        </div>
                    </fieldset>
                </form>
            </div>
            <div class="modal-footer">
                <button type="button" class="btn btn-primary">Save changes</button>
                <button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
            </div>
        </div>
    </div>
</div>
`;

export class Dialog {
    model: mm.Model;
    readonly datatable: DataTables.Api;
    readonly dialog_element: HTMLElement;

    constructor(model: mm.Model) {
        document.body.appendChild(template.content.cloneNode(true));
        this.dialog_element = <HTMLElement>document.body.getElementsByClassName("dialog-reaction-bulkadd")[0]!;
        this.model = model;
    }
}
