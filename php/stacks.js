//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

function toggle_div(id, root_path, page_state_form) {
    var tr_obj = document.getElementById(id); 
    var caret  = document.getElementById(id + "_img");

    if (tr_obj.style.display == "none") {
	tr_obj.style.display = "";
	caret.src            = root_path + "/caret-d.png"; 

	update_page_state(page_state_form, id, 1);
    }
    else {
	tr_obj.style.display = "none";
	caret.src            = root_path + "/caret-u.png";

	update_page_state(page_state_form, id, 0);
    }
}

function toggle_aln_tr(id, root_path, url) {
    var tr_obj = document.getElementById(id); 
    var caret  = document.getElementById(id + "_img");

    if (tr_obj.style.display == "none") {
	tr_obj.style.display = "";
	caret.src            = root_path + "/caret-d.png"; 

	// Check if the alignment has been loaded, if not, load
        // it in the iframe.
        var iframe_obj = document.getElementById(id + '_iframe');
	iframe_obj.src = url;
    }
    else {
	tr_obj.style.display = "none";
	caret.src            = root_path + "/caret-u.png"; 
    }
}

function toggle_sel(id) {
    var div_obj = document.getElementById(id); 

    if (div_obj.style.display == "none") {
        div_obj.style.display = "";
    }
    else {
        div_obj.style.display = "none";
    }
}

function set_operation(id, operation) {

    if (operation == "reset") {
        var verify = confirm("Are you sure you want to delete all corrections for this marker?");
        if (!verify) return;
    }

    var form_obj = document.getElementById(id); 
    form_obj.op.value = operation;

    form_obj.submit();
}

function toggle_cb(form_id, value) {
    var form_obj = document.getElementById(form_id); 

    for(i = 0; i < form_obj.elements.length; i++)
        if (form_obj.elements[i].value == value) 
            if (form_obj.elements[i].checked == true)
                form_obj.elements[i].checked = false;
            else
                form_obj.elements[i].checked = true;
}
