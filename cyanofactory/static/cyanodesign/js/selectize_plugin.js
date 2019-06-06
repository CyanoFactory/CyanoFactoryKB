/**
 * Plugin: "click_edit" (selectize.js)
 * Copyright (c) 2013 Kris Salvador & contributors
 *
 * -- Description
 * - Given an input selectize, users now have the ability to
 * - click on a given item and edit them within the input.
 *
 * example gif: http://g.recordit.co/q391ZyVBVS.gif
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at:
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under
 * the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF
 * ANY KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 *
 * @author Kris Salvador <krissalvador27@gmail.com>
 */

Selectize.define('add_option_hook', function(options) {
	var self = this;

	const registerOptionOrig = this.registerOption;

	const hash_key = function(value) {
		if (typeof value === 'undefined' || value === null) return null;
		if (typeof value === 'boolean') return value ? '1' : '0';
		return value + '';
	};

	this.registerOption = function(data) {
	    var key = hash_key(data[self.settings.valueField]);
	    if (typeof key === 'undefined' || key === null) return false;
	    data.$order = data.$order || ++this.order;
	    this.options[key] = data;

	    /*if (this.options.hasOwnProperty(key)) {
	        return false;
        }*/

	    return key;
    };
});
