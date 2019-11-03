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
    readonly history: HistoryEntry[];
    private readonly app: AppManager;

    constructor(history: HistoryEntry[], app: app.AppManager) {
        this.history = history;
        this.app = app;
    }

    push(entry: any, group: HistoryGroup | null = null) {
        let hentry = new HistoryEntry();

        hentry.type = entry["type"];
        hentry.op = entry["op"];
        hentry.id = entry["id"];
        hentry.object = entry["object"];
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

        if (op == "add") {
            return `${id} (${obj["name"]})`;
        } else if (op == "edit") {
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

            /*return "Edit reaction " + name;
            if (name != obj["name"]) {
                return "Renamed reaction " + name + " to " + obj["name"];
            }*/
        } else if (op == "delete") {
            return `${id} (${obj["name"]})`;
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
        } else if (op == "delete") {
            if (!fn.has("id", id)) {
                throw new Error(entry.type + " " + entry.id + " is not in the model!");
            }

            fn.remove("id", id);
        }
    }

    undo(idx: number): void {
        for (let i = idx; i < this.history.length; ++i) {
            this.history[i].undo = false;
        }

        this.app.model = new mm.Model();
        this.app.model.fromJson(this.app.old_models[0]);

        for (let i = 0; i < this.history.length; ++i) {
            if (!this.history[i].undo) {
                break;
            }

            this.apply(i);
        }

        this.refresh();
    }

    redo(idx: number): void {
        for (let i = 0; i <= idx && i < this.history.length; ++i) {
            if (!this.history[i].undo) {
                this.history[i].undo = true;
                this.apply(i);
            }
        }

        this.refresh();
    }

    refresh(): void {
        this.app.reaction_page.init();
        this.app.metabolite_page.init();
        this.app.settings_page.init();
        this.app.compartment_page.init();

        this.app.history_page.refresh();
    }
}
