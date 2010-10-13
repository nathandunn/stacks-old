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

