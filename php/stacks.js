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
