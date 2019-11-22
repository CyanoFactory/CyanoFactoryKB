import * as app from "./app"
import * as mm from "./metabolic_model";
import * as $ from "jquery";
import {AppManager} from "./app";

export class HistoryEntry {
    type: string;
    op: string;
    id: string;
    object: any;
    undo: boolean = true;
    group: HistoryGroup;

    constructor() {}
}

export class HistoryGroup {
    id: number;
    summary: string;
    date: string;
}

export class HistoryManager {
    history: HistoryEntry[];
    private readonly app: AppManager;
    current_id: number = -1;

    constructor(app: app.AppManager) {
        this.history = [];
        this.app = app;
    }

    push(entry: any, group: HistoryGroup | null = null) {
        let hentry = new HistoryEntry();

        hentry.type = entry["type"];
        hentry.op = entry["op"];
        hentry.id = entry["id"];
        hentry.object = entry["object"];
        if (group == null) {
            group = new HistoryGroup();
            group.id = this.current_id;
            group.date = new Date().toDateString();
            group.summary = "Unsaved changes";
        }
        hentry.group = group;

        let i = 0;
        while (i < this.history.length) {
            if (!this.history[i].undo) {
                this.history.splice(i, 1);
            } else {
                ++i;
            }
         }

        this.history.push(hentry);

        this.app.history_page.refresh();
    }

    op(idx: number): string {
        const op: string = this.history[idx]["op"];
        const first: string = op[0].toUpperCase();
        return first + op.substr(1) + " " + this.history[idx]["type"];
    }

    source(idx: number): string {
        return this.history[idx]["id"];
    }

    parse(idx: number): string {
        const op: string = this.history[idx]["op"];
        const typ: string = this.history[idx]["type"];
        const id: string = this.history[idx]["id"];
        const obj: any = this.history[idx]["object"];

        if (op == "add" || op == "edit") {
            let s = "";

            let first: boolean = true;
            for (const item in obj) {
                if (!first) {
                    s += ", ";
                }
                first = false;

                s += item + ": " + obj[item];
            }
            return s;
        } else if (op == "delete") {
            return `${id}`;
        }
    }

    apply(idx: number): void {
        const entry: HistoryEntry = this.history[idx];

        const op: string = entry.op;
        if (op != "add" && op != "edit" && op != "delete") {
            throw new Error(entry.op + " is not a valid operation!");
        }

        const type: string = entry.type;

        if (type != "reaction" &&
            type != "metabolite" &&
            type != "compartment" &&
            type != "objective") {
            throw new Error(entry.type + " is not a valid target type!");
        }

        const fn: mm.Internal.LstOp<mm.Model> = this.app.model[entry.type];
        const id = entry.id;

        if (op == "add") {
            if (fn.has("id", id)) {
                throw new Error(entry.type + " " + entry.id + " is already in the model!");
            }

            let obj: any = fn.create();
            const eobj: any = entry.object;

            if (id != eobj.id) {
                throw new Error(entry.type + " " + entry.id + " ID mismatch!")
            }

            for (const key in eobj) {
                if (eobj.hasOwnProperty(key)) {
                    if (obj.hasOwnProperty(key)) {
                        obj[key] = eobj[key];
                    }
                }
            }

            obj.fixup();

            fn.add(obj);
        } else if (op == "edit") {
            let obj: any = fn.get("id", id);

            if (obj == null) {
                throw new Error(entry.type + " " + entry.id + " is not in the model!");
            }

            const eobj: any = entry.object;

            if (id != eobj.id && fn.has("id", eobj.id)) {
                throw new Error(entry.type + " " + eobj.id + " is already in the model!");
            }

            for (const key in eobj) {
                if (eobj.hasOwnProperty(key)) {
                    if (obj.hasOwnProperty(key)) {
                        obj[key] = eobj[key];
                    }
                }
            }

            obj.fixup();
        } else if (op == "delete") {
            if (!fn.has("id", id)) {
                throw new Error(entry.type + " " + entry.id + " is not in the model!");
            }

            fn.remove("id", id);
        }
    }

    undo(idx: number): void {
        const undo_func = () => {
            for (let i = 0; i < this.history.length; ++i) {
                if (!this.history[i].undo) {
                    break;
                }

                if (this.current_id == this.history[i].group.id) {
                    console.log("Undo: " + this.history[i]);
                    this.apply(i);
                }
            }

            this.refresh();
        };

        for (let i = idx; i < this.history.length; ++i) {
            this.history[i].undo = false;
        }

        if (this.history[idx].group.id != this.current_id) {
            this.current_id = this.history[idx].group.id;

            // Is at a different savepoint, fetch that model and refresh
            this.app.request_handler.beginRequest();
            $.ajax({
                url: this.app.urls.get_reactions,
                context: document.body,
                data: {
                    revision: this.current_id == -1 ? "" : this.current_id
                }
            }).done((model: any) => {
                this.app.model = new mm.Model();
                this.app.model.fromJson(model);
                this.app.old_model = model;

                undo_func();

                this.app.request_handler.endRequest();
            });
        } else {
            this.app.model = new mm.Model();
            this.app.model.fromJson(this.app.old_model);

            undo_func();
        }
    }

    redo(idx: number): void {
        const redo_func = () => {
            for (let i = 0; i <= idx && i < this.history.length; ++i) {
                if (!this.history[i].undo) {
                    this.history[i].undo = true;

                    if (this.current_id == this.history[i].group.id) {
                        console.log("Redo: " + this.history[i]);
                        this.apply(i);
                    }
                }
            }

            this.refresh();
        };

        if (this.history[idx].group.id != this.current_id) {
            this.current_id = this.history[idx].group.id;

            // Is at a different savepoint, fetch that model and refresh
            this.app.request_handler.beginRequest();
            $.ajax({
                url: this.app.urls.get_reactions,
                context: document.body,
                data: {
                    revision: this.current_id == -1 ? "" : this.current_id
                }
            }).done((model: any) => {
                this.app.model = new mm.Model();
                this.app.model.fromJson(model);
                this.app.old_model = model;

                redo_func();
                this.app.request_handler.endRequest();
            });
        } else {
            redo_func();
        }
    }

    current(): HistoryEntry[] {
        let history = [];
        for (const hist of this.history) {
            if (hist.group.id == this.current_id && hist.undo) {
                history.push(hist);
            }
        }
        return history;
    }

    clear(): void {
        this.current_id = -1;
        this.history = [];
    }

    refresh(): void {
        this.app.reaction_page.init();
        this.app.metabolite_page.init();
        this.app.settings_page.init();
        this.app.compartment_page.init();

        this.app.history_page.refresh();

        this.app.simulation_page.solve();
    }
}
