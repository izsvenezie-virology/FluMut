from PySide6.QtWidgets import QMenu

from flumut_db_editor.database_operations import (
    ForeignKeyViolationError,
)
from flumut_db_editor.gui.dialogs import (
    ConfirmationDialog,
    ErrorDialog,
    ForeignKeyViolationDialog,
    SuccessNotification,
)


class CrudMixin:
    """Mixin to add CRUD operations and context menus to tabs."""

    def show_context_menu(self, pos):
        """Display context menu with edit/delete options."""
        menu = QMenu(self)

        edit_action = menu.addAction('Edit')
        delete_action = menu.addAction('Delete')
        menu.addSeparator()
        refresh_action = menu.addAction('Refresh')

        action = menu.exec(self.mapToGlobal(pos))

        if action == edit_action:
            self.handle_edit()
        elif action == delete_action:
            self.handle_delete()
        elif action == refresh_action:
            self.on_refresh()

    def handle_edit(self):
        """Handle edit action - to be overridden by subclasses."""
        raise NotImplementedError

    def handle_delete(self):
        """Handle delete action - to be overridden by subclasses."""
        raise NotImplementedError

    def on_refresh(self):
        """Refresh data from database."""
        if hasattr(self, 'load_data'):
            self.load_data()


class TableCrudMixin(CrudMixin):
    """CRUD mixin for table-based tabs (BaseTableTab)."""

    def get_selected_item(self):
        """Get the currently selected row in the table."""
        if not hasattr(self, 'table'):
            return None
        current_row = self.table.currentRow()
        if current_row < 0:
            ErrorDialog.show_error(self, 'No Selection', 'Please select an item first.')
            return None
        return current_row

    def handle_edit(self):
        """Edit the selected item - to be overridden."""
        raise NotImplementedError('Subclass must implement handle_edit() with specific model')

    def handle_delete(self):
        """Delete the selected item - to be overridden."""
        raise NotImplementedError('Subclass must implement handle_delete() with specific model')


class TreeCrudMixin(CrudMixin):
    """CRUD mixin for tree-based tabs (HierarchicalTab)."""

    def get_selected_item(self):
        """Get the currently selected item in the tree."""
        if not hasattr(self, 'tree'):
            return None
        current_item = self.tree.currentItem()
        if current_item is None:
            ErrorDialog.show_error(self, 'No Selection', 'Please select an item first.')
            return None
        return current_item

    def get_item_data(self, item):
        """Extract item type and ID from tree item."""
        data = item.data(0, 0x0100)
        if data is None:
            return None, None
        return data

    def handle_edit(self):
        """Edit the selected item - to be overridden."""
        raise NotImplementedError('Subclass must implement handle_edit() with specific model')

    def handle_delete(self):
        """Delete the selected item - to be overridden."""
        raise NotImplementedError('Subclass must implement handle_delete() with specific model')


def delete_with_confirmation(parent, item_type: str, item_id: int, delete_func):
    """
    Common delete workflow: validate → confirm → delete.

    Args:
        parent: Parent widget
        item_type: Type of item (e.g., 'Segment', 'Mutation')
        item_id: ID of the item
        delete_func: Function that performs deletion (must raise ForeignKeyViolationError)
    """
    try:
        # Attempt delete (will raise ForeignKeyViolationError if blocked)
        # But first confirm with user
        if ConfirmationDialog.ask(
            parent,
            f'Delete {item_type}?',
            f'Are you sure you want to delete this {item_type}?',
            'This action cannot be undone.',
        ):
            delete_func(item_id)
            SuccessNotification.show_success(parent, f'{item_type} deleted successfully.')
            return True

    except ForeignKeyViolationError as e:
        ForeignKeyViolationDialog.show_violation(parent, e.item_type, e.item_name, e.violations)

    except Exception as e:
        ErrorDialog.show_error(
            parent,
            'Delete Failed',
            f'Failed to delete {item_type}.',
            str(e),
        )

    return False
